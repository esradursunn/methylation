# methylation
### Analysis of 450k data using minfi

### Installation and Dependencies
rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")

library(minfi)
library(minfiData)
library(sva)

### Reading data
baseDir <- ("/home/edursun/Projects/DKFZ/data/")
list.files(baseDir)
list.files(file.path(baseDir, "5727920038"))
targets <- read.metharray.sheet(baseDir)
targets
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
phenoData[,1:6]
manifest <- getManifest(RGSet)
manifest
head(getProbeInfo(manifest))

### MethylSet and RatioSet
MSet <- preprocessRaw(RGSet)
MSet
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet
BetaValue <- getBeta(RSet)
MValue <- getM(RSet)
CNvalue <- getCN(RSet)

### Save tables
library(writexl)
BetaValue_edit <- cbind(rownames(BetaValue), data.frame(BetaValue, row.names=NULL))
colnames(BetaValue_edit)[colnames(BetaValue_edit) == "rownames(BetaValue)"] <- "Probe"
write_xlsx(BetaValue_edit, path = "/home/edursun/Projects/DKFZ/minfi_results/tables/BetaValue.xlsx", col_names=T)

MValue_edit <- cbind(rownames(MValue), data.frame(MValue, row.names=NULL))
colnames(MValue_edit)[colnames(MValue_edit) == "rownames(MValue)"] <- "Probe"
write_xlsx(MValue_edit, path = "/home/edursun/Projects/DKFZ/minfi_results/tables/MValue.xlsx", col_names=T)

CNvalue_edit <- cbind(rownames(CNvalue), data.frame(CNvalue, row.names=NULL))
colnames(CNvalue_edit)[colnames(CNvalue_edit) == "rownames(CNvalue)"] <- "Probe"
write_xlsx(CNvalue_edit, path = "/home/edursun/Projects/DKFZ/minfi_results/tables/CNvalue.xlsx", col_names=T)

### GenomicRatioSet
GRset <- mapToGenome(RSet)
GRset
beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)
gr <- granges(GRset)
head(gr, n= 3)
annotation <- getAnnotation(GRset)
names(annotation)

### Quality control
## QC plot
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])
qc <- getQC(MSet)
head(qc)
plotQC(qc)
densityPlot(MSet, sampGroups = phenoData$Sample_Group)
densityBeanPlot(MSet, sampGroups = phenoData$Sample_Group)

### Raw samples PCA
PCA_raw <- prcomp(t(getMeth(MSet)), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Cell_type = phenoData$Sample_Group)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) + 
  scale_color_brewer(palette = "Paired")

### To create a shinyMethylSet object
library(shinyMethyl)
myShinyMethylSet <- shinySummarize(RGSet)
runShinyMethyl(myShinyMethylSet)

###  Control probes plot
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
qcReport(RGSet, pdf= "qcReport.pdf")

### Sex prediction
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)
plotSex(getSex(GRset, cutoff = -2))

### Preprocessing and normalization
GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)

### After normalization PCA

PCA_norm <- prcomp(t(getBeta(GRset.quantile)), scale. = FALSE)
percentVar <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG_NORM <- data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2],
                          Cell_type = phenoData$Sample_Group)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) +
  ggtitle("PCA plot of the log-transformed normalized expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) + 
  scale_color_brewer(palette = "Paired")

### Genetic variants and cell type composition
# SNPs
snps <- getSnpInfo(GRset)
head(snps,10)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)

### Cell type composition
library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(RGSet)
cellCounts_edit <- cbind(rownames(cellCounts), data.frame(cellCounts, row.names=NULL))
colnames(cellCounts_edit)[colnames(cellCounts_edit) == "rownames(cellCounts)"] <- "Probe"
write_xlsx(cellCounts_edit, path = "/home/edursun/Projects/DKFZ/minfi_results/tables/cellCounts.xlsx", col_names=T)

### Cell type composition PCA

exp_raw <- log2(cellCounts)
PCA_raw <- prcomp(t(cellCounts), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Cell_type = rownames(dataGG))

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) + 
  ggtitle("PCA plot of the estimated cell type") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) + 
  scale_color_brewer(palette = "Paired")

### Cell type estimation heatmap

annotation_for_heatmap <-
  data.frame(Cell_type = Sample_group)
library(RColorBrewer)
ann_colors <- list(group = brewer.pal(10, name = "Paired"))
names(ann_colors$group) <- unique(Sample_group)
rownames(annotation_for_heatmap) <- colnames(t(cellCounts))
pheatmap(t(cellCounts),
         annotation_col = annotation_for_heatmap,
         annotation_colors = ann_colors,
         scale = "column",
         legend = TRUE,
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         clustering_distance_rows = "manhattan",
         clustering_method = "complete",
         main = "", fontsize_col = 10)


### Identifying DMRs and DMPs
#dmpFinder: to find differentially methylated positions (DMPs)
beta <- getBeta(GRset.quantile)
Sample_group  <- pData(GRset.quantile)$Sample_Group
dmp <- dmpFinder(beta, pheno = Sample_group, type = "continuous")
head(dmp)

### Most variable 100 methylated probes selection
library(dplyr)
dmp_edit <- cbind(rownames(dmp), data.frame(dmp, row.names=NULL))
colnames(dmp_edit)[colnames(dmp_edit) == "rownames(dmp)"] <- "Probe"
dmp_ordered <- arrange(dmp_edit, qval)[1:100,]
GRset.quantile_edit <- getBeta(GRset.quantile)
GRset.quantile_edit_n <- cbind(rownames(GRset.quantile_edit), data.frame(GRset.quantile_edit, row.names=NULL))
colnames(GRset.quantile_edit_n)[colnames(GRset.quantile_edit_n) == "rownames(GRset.quantile_edit)"] <- "Probe"
hundred_genes_raw <- merge(GRset.quantile_edit_n, dmp_ordered, by="Probe")

### dmp Heatmap
hundred_genes_raw.m <- data.frame(hundred_genes_raw, row.names = 1)
ann_col <- data.frame(group = Sample_group)
rownames(ann_col) <- colnames(hundred_genes_raw.m)
Data <- subset(hundred_genes_raw.m, select = -c(intercept, beta, t, pval, qval))
annotation_for_heatmap <-
  data.frame(Cell_type = Sample_group)
library(RColorBrewer)
ann_colors <- list(group = brewer.pal(10, "Paired"))
names(ann_colors$group) <- unique(Sample_group)
rownames(annotation_for_heatmap) <- colnames(Data)
pheatmap(t(Data),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         scale = "row",
         legend = TRUE,
         cluster_cols = T,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = T,
         clustering_distance_rows = "manhattan",
         clustering_method = "complete",
         main = "", fontsize_col = 10)

### bumphunter: to find differentially methylated regions (DMRs)
pheno <- pData(GRset.quantile)$Sample_Group
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(GRset.quantile, design = designMatrix, 
                   cutoff = 0.2, B=0, type="Beta")
library(doParallel)
registerDoParallel(cores = 4)
dmrs <- bumphunter(GRset.quantile, design = designMatrix, 
                   cutoff = 0.2, B=1000, type="Beta")
names(dmrs)
head(dmrs$table, n=3)
data("dmrs_B1000_c02")
head(dmrs$table)

###  Batch effects correction with SVA
library(sva)
mval <- getM(GRset)[1:5000,]
pheno <- pData(GRset)
mod <- model.matrix(~as.factor(status), data=pheno)
mod0 <- model.matrix(~1, data=pheno)
sva.results <- sva(mval, mod, mod0)

###  Other tasks
# A/B compartments prediction
ab <- compartments(grset.quantile, chr="chr14", resolution="100`1000")

# getSnpBeta
snps <- getSnpBeta(RGSet)
head(snps)

#  Out-of-band probes
oob <- getOOB(RGSet)

### Probes in Promoter Region
library(ELMER)
Promoter_probe <- get.feature.probe(TSS, genome = "hg38",
                  met.platform = "450K", TSS.range = list(upstream = 2000, downstream =
                                                            2000), promoter = TRUE, rm.chr = NULL)

### Selecting only probes in promoter region

promoter_names <- data.frame(promoter_probes.edit, colnames("Probe"))
colnames(promoter_names) <- "Probe"
betaValue.n <- data.frame(BetaValue)
betaValue.edit <- cbind(rownames(BetaValue), data.frame(BetaValue, row.names=NULL))
colnames(betaValue.edit)[colnames(betaValue.edit) == "rownames(BetaValue)"] <- "Probe"
promoter_selected <- merge(promoter_names, betaValue.edit, by="Probe")

### PCA for probes in only promoter region
promoter_selected <- data.frame(promoter_selected, row.names = 1)
promoter_selected.na <- na.omit(promoter_selected)
PCA_raw <- prcomp(t(getMeth(MSet)), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Cell_type = promoter_names)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) +
  ggtitle("PCA plot of the probes in promoter region") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) + 
  scale_color_brewer(palette = "Paired")


### BetaValue information for all promoter and non promoter region
library(dplyr)
betaValue_new <- mutate(betaValue.edit, Probe_status =  betaValue.edit$Probe %in% promoter_names$Probe)
Probe_annotation <- betaValue_new$Probe_status
library(plyr)
Probe_annotation <- revalue(as.factor(Probe_annotation), c("TRUE" = "Promoter", "FALSE" = "Non_Promoter"))
Probe_annotation <- data.frame(Probe_annotation)
Probe_annotation.n <- cbind(betaValue_new$Probe, Probe_annotation)
names(Probe_annotation.n)[1] <- "Probe"
names(Probe_annotation.n)[2] <- "Probe_status"


### Heatmap for 100 probes with promoter annotation
names(Probe_annotation)[1] <- "Probe"
dmp_ordered <- arrange(dmp_edit, qval)[1:100,]
hundred_genes_raw_promoter <- merge(hundred_genes_raw, Probe_annotation.n, by="Probe")
hundred_genes_raw_promoter.m <- data.frame(hundred_genes_raw_promoter, row.names = 1)
annotation_for_heatmap <-
  data.frame(Cell_type = hundred_genes_raw_promoter$Probe_status)
ann_col <- data.frame(group = hundred_genes_raw_promoter$Probe_status)
Data_new <- subset(hundred_genes_raw_promoter.m, select = -c(intercept, beta, t, pval, qval, Probe_status))
library(RColorBrewer)
ann_colors <- list(group = c("Promoter" = "green", "Non_promoter" = "darkred"))
names(ann_colors$group) <- unique(hundred_genes_raw_promoter$Probe_status)
rownames(annotation_for_heatmap) <- colnames(t(Data_new))
pheatmap(t(Data_new),
         annotation_col = annotation_for_heatmap,
         annotation_colors = ann_colors,
         scale = "row",
         legend = TRUE,
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F,
         clustering_distance_rows = "manhattan",
         clustering_method = "complete",
         main = "Heatmap for top100 variable CpG's", fontsize_col = 10)


### PCA for 100 Probes with promoter annotation 
PCA_norm <- prcomp((Data_new), scale. = FALSE)
percentVar <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG_norm <- data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2],
                          Cell_type = hundred_genes_raw_promoter$Probe_status)

ggplot(dataGG_norm, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type, shape = Cell_type)) +
  ggtitle("PCA plot of the 100 probes") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_brewer(palette = "Paired")+
  theme().title = »
  
  ### Annotatr
library(annotatr)
library("AnnotationHub")
annots <- c('hg19_genes_promoters')
annots_gr = build_annotations(genome = 'hg19', annotations = annots)
anno <- build_annotations(genome = "hg19", annotations = "hg19_basicgenes")
dm_annotated = annotate_regions(
  regions = gr,
  annotations = annots_gr,
  ignore.strand = TRUE,
  quiet = FALSE)
print(dm_annotated)
df_dm_annotated = data.frame(dm_annotated)
df_dm_annotated_annotr_unique <- unique(df_dm_annotated)
print(head(df_dm_annotated))
annotatr_probes <- dm_annotated@ranges
cpg_annotatr <- annotatr_probes@NAMES
cpg_annotatr.df <- data.frame(cpg_annotatr)

### Selecting cpg in promoter region with annottar
colnames(cpg_annotatr.df) <- "Probe"
betaValue.n <- data.frame(BetaValue)
betaValue.edit <- cbind(rownames(BetaValue), data.frame(BetaValue, row.names=NULL))
colnames(betaValue.edit)[colnames(betaValue.edit) == "rownames(BetaValue)"] <- "Probe"
promoter_selected.annotatr <- merge(cpg_annotatr.df, betaValue.edit, by="Probe")

### PCA for probes in only promoter region
dmp_ordered <- arrange(dmp_edit, qval)
dmp_sorted_genes.annotatr <- merge(promoter_selected.annotatr, dmp_ordered, by="Probe")
dmp_sorted_genes.na.annotatr <- na.omit(dmp_sorted_genes.annotatr)
dmp_sorted_genes.na.annotatr_unique <- unique(dmp_sorted_genes.na.annotatr)
thousand_genes_ordered.annotatr <- arrange(dmp_sorted_genes.na.annotatr_unique, qval)[1:1000,]
thousand_genes_promoter.annotatr <- subset(thousand_genes_ordered.annotatr, select = -c(intercept, beta, t, pval, qval))
thousand_genes_promoter_annotatr.m <- data.frame(thousand_genes_promoter.annotatr, row.names = 1)
PCA_raw_promoter.annotatr <- prcomp(t(thousand_genes_promoter_annotatr.m), scale. = FALSE)
percentVar <- round(100*PCA_raw_promoter.annotatr$sdev^2/sum(PCA_raw_promoter.annotatr$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw_promoter.annotatr$x[,1], PC2 = PCA_raw_promoter.annotatr$x[,2],
                     Cell_type = phenoData$Sample_Group)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) +
  ggtitle("PCA plot of the 1000 most variable probes in promoter region by annotatr") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_brewer(palette = "Paired")

### Homo.sapiens package promoter selection
library(Homo.sapiens)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
proms <- promoters(txdb, upstream=2000,  downstream=200)
dm_annotated_h = annotate_regions(
  regions = gr,
  annotations = proms,
  ignore.strand = TRUE,
  quiet = FALSE)
dm_annotated_h.unique <- unique(dm_annotated_h)
df_dm_annotated_h.unique = data.frame(dm_annotated_h.unique)
print(head(df_dm_annotated_h))
homosapiens_probes <- dm_annotated_h.unique@ranges
cpg_homosapiens <- homosapiens_probes@NAMES
cpg_homosapiens.df <- data.frame(cpg_homosapiens)

### Selecting cpg in promoter region with annottar
colnames(cpg_homosapiens.df) <- "Probe"
betaValue.n <- data.frame(BetaValue)
betaValue.edit <- cbind(rownames(BetaValue), data.frame(BetaValue, row.names=NULL))
colnames(betaValue.edit)[colnames(betaValue.edit) == "rownames(BetaValue)"] <- "Probe"
promoter_selected.homosapiens <- merge(cpg_homosapiens.df, betaValue.edit, by="Probe")

### PCA for probes in only promoter region
dmp_ordered <- arrange(dmp_edit, qval)
dmp_sorted_genes.homosapiens <- merge(promoter_selected.homosapiens, dmp_ordered, by="Probe")
dmp_sorted_genes.na.homosapiens <- na.omit(dmp_sorted_genes.homosapiens)
dmp_sorted_genes.na.homosapiens_unique <- unique(dmp_sorted_genes.na.homosapiens)
thousand_genes_ordered.homosapiens <- arrange(dmp_sorted_genes.na.homosapiens_unique, qval)[1:1000,]
thousand_genes_promoter.homosapiens <- subset(thousand_genes_ordered.homosapiens, select = -c(intercept, beta, t, pval, qval))
thousand_genes_promoter_homosapiens.m <- data.frame(thousand_genes_promoter.homosapiens, row.names = 1)
PCA_raw_promoter.homosapiens <- prcomp(t(thousand_genes_promoter_homosapiens.m), scale. = FALSE)
percentVar <- round(100*PCA_raw_promoter.homosapiens$sdev^2/sum(PCA_raw_promoter.homosapiens$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw_promoter.homosapiens$x[,1], PC2 = PCA_raw_promoter.homosapiens$x[,2],
                     Cell_type = phenoData$Sample_Group)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell_type)) +
  ggtitle("PCA plot of the 1000 most variable probes in promoter region by Homo.sapiens") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_brewer(palette = "Paired")

