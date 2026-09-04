library("DESeq2")

args <- commandArgs(TRUE);
args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/kallisto_output/kallistoData.Rdata", "WT", "cmt3")

load(args[1])
colD <- subset(colD, subset = condition %in% args[c(2,3)]) # subsetting to get only control & treatment conditions
countsToUse <- countsToUse[,colnames(countsToUse) %in% colD$sample] # subsetting matrix to get only control & treatment conditions
metadata <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="")

dds = DESeqDataSetFromMatrix(countsToUse,colData=colD,design=~condition)

TEGs <- metadata[metadata$Type=="transposable_element_gene",]

# adding meta data to the dataframe
mcols(dds) <- DataFrame(mcols(dds), metadata)
mcols(dds)
# pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# setting the reference treatment
dds$condition <- relevel(dds$condition, ref = "WT")
# differential analysis
dds <- DESeq(dds)

w2 <- results(dds, contrast=c("condition","w2","WT"))
suvh456 <- results(dds, contrast=c("condition","suvh456","WT"))
w2_suvh456 <- results(dds, contrast=c("condition","w2_suvh456","WT"))
cmt3 <- results(dds, contrast=c("condition","cmt3","WT"))
w2_cmt3 <- results(dds, contrast=c("condition","w2_cmt3","WT"))
h1 <- results(dds, contrast=c("condition","h1","WT"))
w2_h1_F5 <- results(dds, contrast=c("condition","w2_h1_F5","WT"))
w2_h1_H6 <- results(dds, contrast=c("condition","w2_h1_H6","WT"))

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
upTEGs <- function(res){
  TEGs <- subset(metadata, Type=="transposable_element_gene")
  nrow(res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,])
  }
samples <- list(w2, suvh456, w2_suvh456, cmt3, w2_cmt3, h1, w2_h1_F5, w2_h1_H6)
names(samples) <- c("w2", "suvh456", "w2_suvh456", "cmt3", "w2_cmt3", "h1", "w2_h1_F5", "w2_h1_H6")
lapply(samples, FUN=upTEGs)
