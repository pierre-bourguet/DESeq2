#### DESeq2 ON KALLISTO COUNTS - ms1xerj ovules

library("DESeq2")
library("rtracklayer")
library("gProfileR")
library("GeneOverlap")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")

#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 1. Load data from pipeline
setwd("/groups/berger/user/michael.borg/RNA-seq/demeter/results_myRuns/kallisto_output/")
load("kallistoData.Rdata")
##The object **colD** is the colData for DEseq (colD$condition contains the condition for each sample)
##The **tpm** object contains the TPM values for each samples
##The **countsToUse** is the matrix with counts for DEseq

#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 2. Make sure the samples are in the same order in colD and tpm (they should be for the pipeline)
if (! all.equal(colD$sample,colnames(tpm))==T)
  stop("names not matching")


#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 3. Calculate tpm means per condition
tpm_mean = matrix(NA, ncol = length(unique(colD$condition)), nrow= nrow(tpm), dimnames = list(rownames(tpm), unique(colD$condition)))
for (c in unique(colD$condition))
  tpm_mean[,c] = rowMeans(tpm[,colD$condition==c])


#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 4. Test that the calculations were okay (using just one hand picked example)
if (! all.equal(tpm_mean[100,1],mean(tpm[100,1:3]))==T)
  stop("calculations not correct")


#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 5. DESeq2 analysis and contrasts
countsToUse <- subset(countsToUse, rowMax(countsToUse) >= 10)
dds <- DESeqDataSetFromMatrix(countsToUse, colData=colD, design=~condition)
filter.dds <- DESeqDataSetFromMatrix(countsToUse, colData = colData(dds), design=~condition)
run.dds <- DESeq(filter.dds)
results <- results(run.dds)


#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 6. PCA plot and cor.matrix
# Regularized log transformation for different analysis (clustering, heatmaps, etc)
rld <- rlog(run.dds, blind=FALSE)
vst <- vst(run.dds, blind=FALSE)

pdf("/groups/berger/user/michael.borg/RNA-seq/demeter/results_myRuns/pca_rld.pdf", width = 5, height = 5)
plotPCA(rld, intgroup=c("condition"))
dev.off()

pdf("/groups/berger/user/michael.borg/RNA-seq/demeter/results_myRuns/pca_vst.pdf", width = 5, height = 5)
plotPCA(vst, intgroup=c("condition"))
dev.off()

cor.matrix <- cor(tpm, method = "spearman")
pheatmap(cor.matrix, 
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = TRUE,
         clustering_distance_cols = "correlation",
         show_rownames = TRUE,
         scale = "none",
         border_color = NA, 
         main = "",
         cutree_rows = 2,
         cutree_cols = 2,
         #color = colours,
         width = 7,
         height = 6)

#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- 7. Perform the DESEQ comparisons
## import annotations
setwd("/groups/berger/user/michael.borg/ATAC-seq/gene_lists/Araport/")
ara11 <- import.gff3("Araport11_GFF3_genes_transposons.201606.gff")
score(ara11[is.na(score(ara11))]) <- 0
seqlevels(ara11) <- c("1","2","3","4","5","Pt","Mt")
ara11.df <- as.data.frame(subset(ara11, locus_type == "protein_coding"))

########### --------------
tpm_mean <- as.data.frame(tpm_mean)
tpm_mean$ID <- rownames(tpm_mean)

setwd("/groups/berger/user/michael.borg/RNA-seq/demeter/results_myRuns/")
contrasts <- read.table("contrasts.txt", header = T, sep=",")
DEG.list <- list()
for (i in 1:nrow(contrasts)){
  co <- c("condition", colnames(contrasts)[c(which(contrasts[i,]==1), which(contrasts[i,]==-1))])
  x <- as.data.frame(results(run.dds, contrast = co))  ## use <lfcShrink> instead to give the MLE log2FC
  x[,7] <- rownames(x)
  colnames(x)[7] <- "ID"
  z <- merge(x, ara11.df, by = "ID")
  z <- merge(z, tpm_mean, by = "ID")
  DEG.list[[i]] <- z[,c(1:7,34,33,17:21,30)]
  y <- paste(co[2], "-vs-", co[3], sep = "")
  names(DEG.list)[[i]] <- y
}
rm(x,y,z,i,contrasts)

DEGs <- DEG.list[[1]]

setwd("/groups/berger/user/michael.borg/RNA-seq/demeter/results_myRuns/")
write.table(DEGs, file=paste(names(DEG.list)[[1]], "txt", sep="."), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
