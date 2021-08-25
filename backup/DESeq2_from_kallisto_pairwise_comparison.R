#!/usr/bin/env Rscript 
library("DESeq2")
args <- commandArgs(TRUE);
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/kallisto_output/kallistoData.Rdata", "WT", "cmt3", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/", 3)

setwd(args[4])
load(args[1])
colD_subset <- subset(colD, subset = condition %in% args[c(2,3)]) # subsetting to get only control & treatment conditions
countsToUse_subset <- countsToUse[,colnames(countsToUse) %in% colD_subset$sample] # subsetting matrix to get only control & treatment conditions
counts_RPM <- read.delim("kallisto_output/kallisto_tpms.tab", head=T, sep=" ", quote="", comment.char="")
names(counts_RPM) <- paste(colD$condition, "_rep", c(1:args[5]), sep="")
row.names(counts_RPM) <- gsub(x=row.names(counts_RPM), pattern= "\"", replacement="")
metadata <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="")
row.names(metadata) <- metadata$GeneId
dds = DESeqDataSetFromMatrix(countsToUse_subset,colData=colD_subset,design=~condition)

# adding meta data to the dataframe
# mcols(dds) <- DataFrame(mcols(dds), metadata) # doesn't work since the number of annotations is different in both files, I would have to create a new one, but anyway i didn't use this metadata so far
# pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# setting the reference treatment
dds$condition <- relevel(dds$condition, ref = args[2])
# differential analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", args[3], args[2]))

# DEGs ####
TEGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element_gene")
PCGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="gene")
upTEGs <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
upPCGs <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
downTEGs <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
downPCGs <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
DEGs <- list(upTEGs, upPCGs, downTEGs, downPCGs)
names(DEGs) <- c("upTEGs", "upPCGs", "downTEGs", "downPCGs")
# functions to normalize raw counts to RPM
df <- upTEGs
fill_data_to_df <- function(df, sample_nb) {
  df <- merge(as.data.frame(df), metadata, by="row.names")
  row.names(df) <- df[,1]
  df <- merge(df[,-1], counts_RPM, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  df$Length <- df$End - df$Start + 1
  return(df[,c(8:11,(14+sample_nb),12,7,1:6,14:(14+sample_nb-1),13)])
}
colnames(countsToUse) <- colD$condition
DEGs <- lapply(DEGs, FUN=fill_data_to_df, sample_nb=length(colD$condition))
# export DEG tables
dirname <- paste(args[3], "_vs_", args[2], "_DEGs/", sep="")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
mapply(FUN=write.table, x=DEGs, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)
