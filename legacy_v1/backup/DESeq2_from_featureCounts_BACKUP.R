
# args <- c("Y:/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/"
          , "counts_summary.tsv"
          , "Y:/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/samples.txt"
          , "Y:/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF")

args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/"
          , "counts_summary.tsv"
          , "/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/samples.txt"
          , "/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF")

setwd(args[1])
# importing DFs
cts_summary <- read.delim(args[2], header=T, sep='\t', quote="", dec=".", comment.char="")
cts <- ceiling(cts_summary[,c(9:17)]) ; row.names(cts) <- cts_summary$Geneid # ceiling is to round up
samples <- read.delim(args[3], header=F, sep='\t', quote="", dec=".", comment.char="")[,1]
metadata <- read.delim(args[4], head=T, sep="\t", quote="", comment.char="")
# defining metadata
coldata <- data.frame(
  condition=gsub("_rep[0-9]{1}", "",samples),
  type=rep("single-read", length(samples))
)
row.names(coldata) <- samples

all(rownames(coldata) == colnames(cts)) # IF FALSE YOU ARE IN DEEP SHIT MY MAN, NOT GONNA WORK

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

# adding meta data to the dataframe
mcols(dds) <- DataFrame(mcols(dds), metadata)
mcols(dds)
# pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# setting the reference treatment
dds$condition <- relevel(dds$condition, ref = "Col0")
# differential analysis
dds <- DESeq(dds)
res_ddm1_7th <- results(dds, contrast=c("condition","ddm1_7th","Col0"))
res_ddm1_1st <- results(dds, contrast=c("condition","ddm1_1st","Col0"))
res <- res_ddm1_1st
# plot a single gene (here the one with the lowest pval)
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
# plot PCA
plotPCA(vsd, intgroup=c("condition", "type"))

# upregulated TEs
nrow(res_ddm1_1st[mcols(dds)$Type == "transposable_element_gene" & res_ddm1_1st$log2FoldChange >=1 & !is.na(res_ddm1_1st$padj) & res_ddm1_1st$padj < 0.05,])
nrow(res_ddm1_7th[mcols(dds)$Type == "transposable_element_gene" & res_ddm1_7th$log2FoldChange >=1 & !is.na(res_ddm1_7th$padj) & res_ddm1_7th$padj < 0.05,])
# upregulated TEs
upTEGs <- function(res){
  TEGs <- subset(read.delim("Y:/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element_gene")
  nrow(res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,])
}
samples <- list(w2, suvh456, w2_suvh456, cmt3, w2_cmt3, h1, w2_h1_F5, w2_h1_H6)
names(samples) <- c("w2", "suvh456", "w2_suvh456", "cmt3", "w2_cmt3", "h1", "w2_h1_F5", "w2_h1_H6")
lapply(samples, FUN=upTEGs)