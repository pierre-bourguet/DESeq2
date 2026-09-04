#!/usr/bin/env Rscript 

args <- commandArgs(TRUE);
# args <- c("Y:/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/", "counts_summary.tsv", "Y:/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/samples.txt", "Y:/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF"))
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/", "Col0", "ddm1_1st", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/samples.txt")
setwd(args[1])
# importing DFs
cts_summary <- read.delim("counts_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="")
row.names(cts_summary) <- cts_summary$Geneid
sample_columns <- c(grep(args[2], names(cts_summary)), grep(args[3], names(cts_summary))) # retrieve columns that contain sample names
samples <- read.delim(args[4], header=F, sep='\t', quote="", dec=".", comment.char="")[,1]
subsamples <- samples[c(grep(args[2], samples), grep(args[3], samples))] # retrieve columns that contain sample names
metadata <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_lncRNA.SAF", head=T, sep="\t", quote="", comment.char=""),  subset = !Chr %in% c("chrC", "chrM"))

DESeq2_function <- function(x) { # x should be cts_summary.tsv file (summary of raw counts)
  x <- cts_summary
  cts <- ceiling(x[,sample_columns]) ; row.names(cts) <- x$Geneid # ceiling is to round up
  cts <- cts[!(substr((row.names(cts)),3,3) == "C" | substr((row.names(cts)),3,3) == "M"),] # removing chloroplastic and mitochondrial genes
  # defining metadata
  coldata <- data.frame(
    condition=gsub("_rep[0-9]{1}", "", subsamples),
    type=rep("single-read", length(subsamples))
  )
  row.names(coldata) <- subsamples
  if (!all(rownames(coldata[1:6,]) == colnames(cts))) {
    print("the names in count_file and the sample table do not match")
  } # IF FALSE YOU ARE IN DEEP SHIT MY MAN, NOT GONNA WORK
  library("DESeq2")
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata[1:6,],
                                design = ~ condition)
  
  # adding meta data to the dataframe
  mcols(dds) <- DataFrame(mcols(dds))
  # pre-filtering, here keeping only
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  # setting the reference treatment
  dds$condition <- relevel(dds$condition, ref = args[2])
  # differential analysis
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition",args[3],args[2]))
  return(res)
}
res <- DESeq2_function(cts_summary)

# DEGs ####
lncRNAs <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_lncRNA.SAF", head=T, sep="\t", quote="", comment.char="")
uplncRNAs <- res[row.names(res) %in% lncRNAs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
downlncRNAs <- res[row.names(res) %in% lncRNAs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
DEGs <- list(uplncRNAs, downlncRNAs)
names(DEGs) <- c("uplncRNAs", "downlncRNAs")
# functions to normalize raw counts to RPM
RPM <- function(prefix, df){ # normalize with number of reads mapped to chr1-5
  x <- read.delim(paste("idxstats/", prefix, "-idxstats.tsv", sep=""), header=T, sep='\t', quote="", dec=".", comment.char="")
  totreads <- sum(subset(x, subset = substr(chr,4,4) %in% c(1:5))$mapped_reads) / 1000000
  df[,which(names(df)==prefix)] <- eval(parse(text=paste("df", "$", prefix, sep=""))) / totreads # eval(parse) converts a character string into an interpretable object AFAIK
  #df[,paste(prefix, "_rpkm", sep="")] <- eval(parse(text=paste("df", "$", prefix, sep=""))) / totreads # eval(parse) converts a character string into an interpretable object AFAIK
  #names(df)[which(names(df)==prefix)] <- paste(prefix, "_rpm", sep="")
  names(df)[which(names(df)==prefix)] <- paste(prefix, "_rpm", sep="")
  return(df)
}
RPM_to_df <- function(df, sample_nb) {
  df <- merge(as.data.frame(df), cts_summary, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  for (i in names(df)[15:(15+sample_nb-1)]) {
    df <- RPM(i, df) # apply RPM normalization to all columns
  }
  return(df[,c(8:13,7,1:6,15:(15+sample_nb-1),14)])
}
DEGs <- lapply(DEGs, FUN=RPM_to_df, sample_nb=length(samples))
# export DEG tables
dirname <- paste(args[3], "_vs_", args[2], "_DEGs/", sep="")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
mapply(FUN=write.table, x=DEGs, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)

#### DEGs in antisense, filtering out DEGs that are already misregulated in the same direction in sense
cts_summary_AS <- read.delim("counts_AS_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="")
res_AS <- DESeq2_function(cts_summary_AS)
# DEGs with AS ####
uplncRNAs_AS <- res_AS[row.names(res_AS) %in% lncRNAs$GeneId & res_AS$log2FoldChange >=1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(uplncRNAs)),]
downlncRNAs_AS <- res_AS[row.names(res_AS) %in% lncRNAs$GeneId & res_AS$log2FoldChange <=-1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(downlncRNAs)),]
DEGs_AS <- list(uplncRNAs_AS, downlncRNAs_AS)
names(DEGs_AS) <- c("uplncRNAs_AS", "downlncRNAs_AS")
DEGs_AS <- lapply(DEGs_AS, FUN=RPM_to_df, sample_nb=length(samples))
# export DEG tables
dirname <- paste(args[3], "_vs_", args[2], "_DEGs_AS/", sep="")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
mapply(FUN=write.table, x=DEGs_AS, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs_AS), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs_AS, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_DEGs_AS_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)