Tuto DADA2
================

# Partie 1 : Vérification des données

## Téléchargement des packages adéquats

`{bash include=FALSE} sudo apt-get update -y sudo apt-get install -y libglpk-dev  sudo apt-get install -y liblzma-dev libbz2-dev`

`{r include=FALSE} if (!requireNamespace("BiocManager", quietly = TRUE))     install.packages("BiocManager") BiocManager::install("BiocStyle") BiocManager::install("Rhtslib")`

``` {r}
library("knitr")
```

``` {r}
library("BiocStyle")
```

``` {r}
.cran_packages <- c("ggplot2", "gridExtra", "devtools")
install.packages(.cran_packages) 
```

``` {r}
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
BiocManager::install(.bioc_packages)
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

`{bash, echo=FALSE} cd ~ wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip unzip miseqsopdata.zip`

`{r, echo=FALSE} install.packages("gitcreds")`

`{r, echo=FALSE} set.seed(100) miseq_path <- "/home/rstudio/MiSeq_SOP" list.files(miseq_path)`

##Lecture des reads

``` {r}
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
```

## Vérification de la qualité des reads

``` {r}
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
```

### Visualisation de la qualité des reads

`{r plotqualityprofile} fnFs[1:3] plotQualityProfile(fnFs[1:2])`

### Filtration des reads et suppression des séquences dont le score qualité est trop faible

``` {r}
plotQualityProfile(fnRs[1:2])
```

``` {r}
filt_path<- file.path(miseq_path,"filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs<-file.path(filt_path, paste0 (sampleNames, "_F_filt.fastq.gz"))
filtRs<-file.path(filt_path, paste0(sampleNames,"_R_filt.fastq.gz"))
 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

### Conversion des séquences identiques en séquence unique

``` {r}
derepFs <-derepFastq(filtFs, verbose = TRUE)
derepRs <-derepFastq(filtRs, verbose = TRUE)
```

``` {r}
names(derepFs)<-sampleNames
names(derepRs)<-sampleNames
```

### Estimation des erreurs de séquençage

``` {r}
errF<-learnErrors(filtFs, multithread = TRUE)
errR<-learnErrors(filtRs, multithread = TRUE)
```

### Visualisation des erreurs de séquençage

``` {r}
plotErrors(errF)
plotErrors(errR)
```

``` {r}
dadaFs<- dada(derepFs, err = errF, multithread = TRUE)
dadaRs<-dada(derepRs, err=errR,multithread = TRUE)
dadaFs[[1]]
```

## Alignement des reverses et forwards

### Création d’une table de séquences

``` {r}
mergers<-mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll<-makeSequenceTable(mergers[!grepl("Mock",names(mergers))])
table((nchar(getSequences(seqtabAll))))
```

### Suppression des chimères

``` {r}
seqtabNoC<- removeBimeraDenovo(seqtabAll)
```

##Attribution de la taxonomie

### Attribution à l’aide de l’assigneur Silva

``` {bash}
cd ~
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
```

``` {r}
fastaRef<-"/home/rstudio/silva_nr99_v138.1_train_set.fa.gz"
taxTab<-assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread = TRUE)
unname(head(taxTab))
```

###construction d’un abre phylogénétique

`{r Alignement des séquences} seqs<- getSequences(seqtabNoC) names(seqs)<-seqs alignment<-AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=FALSE)`

``` {r}
phangAlign<-phyDat(as(alignment, "matrix"), type="DNA")
dm<-dist.ml(phangAlign)
treeNJ<-NJ(dm)
fit= pml(treeNJ, data=phangAlign)
fitGTR<- update(fit,k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

plot(fitGTR)
```

## Exploitation du package phyloseq

### Attribution des données à un objet phyloseq

``` {r}
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] 
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) 
all(rownames(seqtabAll) %in% samdf$SampleID)
```

``` {r}
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]
```

# Partie 2 : Analyse des données à l’aide u package phyloseq

## Chargement des data dans phyloseq

``` {r}
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) 
```

``` {r}
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

## Filtration des données taxonomiques

``` {r}
rank_names(ps)
```

## Création d’un tableau des nombres de lectures pour chaque embranchement

``` {r}
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

## Suppression des données ambiguës

``` {r}
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

``` {r}
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

## Calcul des prévalences totales et moyennes dans chaque embranchement

`{r echo=TRUE} plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})`

## Exploration des embranchements comprenant Deinococcus thermus dans le phylum Fusobacteria

### Filtration des phyla :

`{r echo=TRUE} filterPhyla = c("Fusobacteria", "Deinococcus-Thermus") ps1 = subset_taxa(ps, !Phylum %in% filterPhyla) ps1`

### Filtrage de la prévalence :

`{r echo=TRUE} prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum")) ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +   # Include a guess for parameter   geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +   facet_wrap(~Phylum) + theme(legend.position="none")`

Calcul de la prévalence, donc la présence dans un grand nombre
d’échantillons de tel phyla.

``` {r}
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

``` {r}
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

\`\`\`{r warning=FALSE, include=FALSE} .cran_packages \<- c(
“shiny”,“miniUI”, “caret”, “pls”, “e1071”, “ggplot2”, “randomForest”,
“dplyr”, “ggrepel”, “nlme”, “devtools”, “reshape2”, “PMA”, “structSSI”,
“ade4”, “ggnetwork”, “intergraph”, “scales”) .github_packages \<-
c(“jfukuyama/phyloseqGraphTest”) .bioc_packages \<- c(“genefilter”,
“impute”)


    ```{r include=FALSE}
    .cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                      "reshape2", "PMA", "structSSI", "ade4",
                      "ggnetwork", "intergraph", "scales")
    .github_packages <- c("jfukuyama/phyloseqGraphTest")
    .bioc_packages <- c("genefilter", "impute")

`{r include=FALSE} .inst <- .cran_packages %in% installed.packages() if (any(!.inst)){   install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/") }`

`{r include=FALSE} .inst <- .github_packages %in% installed.packages() if (any(!.inst)){   devtools::install_github(.github_packages[!.inst]) }`

`{r include=FALSE} .inst <- .bioc_packages %in% installed.packages() if(any(!.inst)){BiocManager::install(.bioc_packages[!.inst]) }`

### Combinaison des caractéristiques des bactéries du même genre

``` {r}
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

### Construction d’un arbre d’agglomération

``` {r}
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

`{r echo=TRUE} multiPlotTitleTextSize = 15 p2tree = plot_tree(ps2, method = "treeonly",                    ladderize = "left",                    title = "Before Agglomeration") +   theme(plot.title = element_text(size = multiPlotTitleTextSize)) p3tree = plot_tree(ps3, method = "treeonly",                    ladderize = "left", title = "By Genus") +   theme(plot.title = element_text(size = multiPlotTitleTextSize)) p4tree = plot_tree(ps4, method = "treeonly",                    ladderize = "left", title = "By Height") +   theme(plot.title = element_text(size = multiPlotTitleTextSize)) grid.arrange(nrow = 1, p2tree, p3tree, p4tree)`

### Création d’un graphique d’abondance relative

\`\`\`{r echo=TRUE} plot_abundance = function(physeq,title = ““, Facet
=”Order”, Color = “Phylum”){ # Arbitrary subset, based on Phylum, for
plotting p1f = subset_taxa(physeq, Phylum %in% c(“Firmicutes”)) mphyseq
= psmelt(p1f) mphyseq \<- subset(mphyseq, Abundance \> 0) ggplot(data =
mphyseq, mapping = aes_string(x = “sex”,y = “Abundance”, color = Color,
fill = Color)) + geom_violin(fill = NA) + geom_point(size = 1, alpha =
0.3, position = position_jitter(width = 0.3)) + facet_wrap(facets =
Facet) + scale_y\_log10()+ theme(legend.position=“none”) }

ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})


    ### Visualisation des valeurs d'abondance avant et après transformation

    ```{r echo=TRUE}
    plotBefore = plot_abundance(ps3,"")
    plotAfter = plot_abundance(ps3ra,"")
    # Combine each plot into one graphic.
    grid.arrange(nrow = 2,  plotBefore, plotAfter)

### Abondances relatives de Lactobacilliales

On remarque que les Lactobacilliales est un ordre avec un profil
d’abondance binomial. Recherche d’une explication taxonomique.

`{r echo=TRUE} psOrd = subset_taxa(ps3ra, Order == "Lactobacillales") plot_abundance(psOrd, Facet = "Genus", Color = NULL)`

### Création d’une variable catégorielle correspondant à l’âge des souris

`{r echo=TRUE} qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")`

Trois groupes distincts sont observables

### Comparaison de la profondeur des lectures brutes et transformées en log

`{r echo=TRUE} qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +   xlab("Logged counts-per-sample")`

### Analyse des coordonnées principales (PCoA) avec la dissemblance de Bray-Curtis sur la distance Unifrac pondérée

``` {r}
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

### Vérification des valeurs aberrantes féminines

``` {r}
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

## Calcul des ordinations avec ces valeurs aberrantes supprimées

``` {r}
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

### Suppression des reads avec moins de 1000 lectures

``` {r}
which(!rowSums(otu_table(ps)) > 1000)
```

``` {r}
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

### Création de la PCoA avec dissimilarité de Bray-Curtis

``` {r}
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

### Analyse des coordonnées principales doubles (DPCoA)

``` {r}
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

``` {r}
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

``` {r}
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

### Création d’une PCA sur rangs

``` {r}
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

### 

``` {r}
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

``` {r}
library(dplyr)
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

### Réalisation de l’ACP avec annotations

``` {r}
library(ade4)
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
```

``` {r}
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
col_scores <- col_scores %>%
  left_join(tax)
```

### Réalisation d’un biplot après transformation de rangs tronqués

``` {r}
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

``` {r}
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

``` {r}
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

## Analyse canonique des correspondances (CCpnA)

``` {r}
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

\`\`\`{r echo=TRUE} library(ggrepel) ps_scores \<-
vegan::scores(ps_ccpna) sites \<-
data.frame(ps_scores*s**i**t**e**s*)*s**i**t**e**s*SampleID \<-
rownames(sites) sites \<- sites %>% left_join(sample_data(ps))

species \<-
data.frame(ps_scores*s**p**e**c**i**e**s*)*s**p**e**c**i**e**s*otu_id
\<- seq_along(colnames(otu_table(ps))) species \<- species %>%
left_join(tax) evals_prop \<- 100 \* ps_ccpna*C**C**A*eig\[1:2\] /
sum(ps_ccpna*C**A*eig) ggplot() + geom_point(data = sites, aes(x = CCA1,
y = CCA2), shape = 2, alpha = 0.5) + geom_point(data = species, aes(x =
CCA1, y = CCA2, col = Order), size = 0.5) + geom_text_repel(data =
species %>% filter(CCA2 \< -2), aes(x = CCA1, y = CCA2, label = otu_id),
size = 1.5, segment.size = 0.1) + facet_grid(. \~ family_relationship) +
guides(col = guide_legend(override.aes = list(size = 3))) + labs(x =
sprintf(“Axis1 \[%s%% variance\]”, round(evals_prop\[1\], 2)), y =
sprintf(“Axis2 \[%s%% variance\]”, round(evals_prop\[2\], 2))) +
scale_color_brewer(palette = “Set2”) +
coord_fixed(sqrt(ps_ccpna*C**C**A*eig\[2\] /
ps_ccpna*C**C**A*eig\[1\])\*0.45 ) + theme(panel.border =
element_rect(color = “#787878”, fill = alpha(“white”, 0)))


    ## Recherche d'un biomarqueur à l'aide du machine learning

    ### Division des données en ensemble d'apprentissages et de tests

    ```{r}
    library(caret)
    sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
    dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
    trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
    inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
    training <- dataMatrix[inTrain,]
    testing <- dataMatrix[-inTrain,]
    plsFit <- train(age ~ ., data = training,
                    method = "pls", preProc = "center")

``` {r}
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

### Utilisation des forêts aléatoires pour déterminer le.s biomarqueur.s

``` {r}
library(randomForest)
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

### Réalisation d’un bi-plot associé à des tracés de proximité projettant les échantillons en fonction de l’âge des souris

``` {r}
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

``` {r}
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

### Création d’un graphique de proximité, basé sur la distance calculée entre échantillons estimée par le package random forest

``` {r}
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

`{r echo=TRUE, message=TRUE} as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])`

Le microorganisme du genre Roseburia, appartenant à la famille des
Lachnospiraceae est le biomarqueur le plus discriminant identifié par le
package random forest

``` {r}
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

## Représentations graphiques des résulats

### Création d’un réseau basé sur une distance de Jaccard à 0,35

``` {r}
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]
```

``` {r}
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

``` {r}
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

``` {r}
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```

``` {r}
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

``` {r}
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```

``` {r}
library("nlme")
library("reshape2")
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
#                                  diversity_means$host_subject_id)
```

``` {r}
library("reshape2")
library("DESeq2")
#New version of DESeq2 needs special levels
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
```
