#!/usr/bin/env Rscript 

args <- commandArgs(TRUE);
# args[1] is path to counts_summary.tsv files
# args[2] is reference sample
# args[3] is treatment sample
# args[4] is path to sample table
# args[5] is paired-end or single-strand
# args[6] is number of replicates

# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3/w2_h1_suvh456_cmt3_quantseq/", "WT", "w2_suvh456", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_quantseq.samples", "single-strand", 3)
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3/w2_h1_suvh456_cmt3_polyA//", "WT", "w2_suvh456", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_polyA.samples", "paired-end", 3)
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/adcp1_Zhao_2019/", "WT", "adcp1", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/adcp1.samples", "paired-end", 2)


# importing DFs
setwd(args[1])
nb_replicate <- as.numeric(args[6])
TEGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element_gene")
TEs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element")
PCGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="gene")
cts_summary <- read.delim("counts_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="")
row.names(cts_summary) <- cts_summary$Geneid
cts_summary <- cts_summary[!(substr((row.names(cts_summary)),3,3) == "C" | substr((row.names(cts_summary)),3,3) == "M"),] # removing chloroplastic and mitochondrial genes
samples <- read.delim(args[4], header=T, sep=',', quote="", dec=".", comment.char="")[,c(2,3)]
samples <- samples[order(samples[,1]),]
subsamples <- subset(samples, samples$condition %in% c(args[2], args[3]))
subsample_columns <- which(names(cts_summary) %in% paste(subsamples$condition, subsamples$sample, sep="_")) # retrieve columns that contain subsample names
sample_columns <- which(names(cts_summary) %in% paste(samples$condition, samples$sample, sep="_")) # retrieve columns that contain sample names
# metadata <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""),  subset = !Chr %in% c("chrC", "chrM")) # not using this so far

# import PCG exon size to perform TPM normalization on cDNA length (for PCG only, TEs are normalized on full length)
exon_size <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_gtf_exonic_genes_sizes.tsv", head=F, sep="\t", quote="", comment.char="", col.names=c("Geneid", "Length"))
tmp <- merge(cts_summary, exon_size, by="Geneid", sort=F) # input files should already be in the same order
if (sum(droplevels(cts_summary$Geneid[cts_summary$Geneid %in% tmp$Geneid]) != droplevels(tmp$Geneid))>0) {
  print("re-ordering of dataframes went wrong")
} # check that ordering is correct
cts_summary_TPM <- cts_summary
cts_summary_TPM$Length[cts_summary_TPM$Geneid %in% tmp$Geneid] <- tmp$Length.y # replaces length by exon length

# TPM normalization
TPM <- function(x) { # TPM normalization function
  return( (x / (cts_summary_TPM$Length / 1000)) / (sum(x / (cts_summary_TPM$Length / 1000)) / 1000000) )
}
cts_summary_TPM[,sample_columns] <- apply(X=cts_summary_TPM[,sample_columns], MARGIN = 2, FUN = TPM) # normalizes all columns
# creates dataframe with average values among 3 replicates
rep_mean <- function(x) {
  y <- as.data.frame(apply(array(as.matrix(x[sample_columns]), c(nrow(x), 3, length(unique(samples$condition)))), 3, rowMeans, na.rm = TRUE)) # calculate the average for every 3 columns (assuming 3 rep) and returns a matrix that we convert to DF
  names(y) <- unique(samples$condition)
  return(y)
}
cts_summary_TPM_mean <- cbind(cts_summary_TPM[which(!names(cts_summary) %in% names(cts_summary)[sample_columns])]
                              , rep_mean(cts_summary_TPM))

# DESeq2 function
suppressPackageStartupMessages(library("DESeq2"))
DESeq2_function <- function(x) { # x should be cts_summary.tsv file (summary of raw counts)
  cts <- ceiling(x[,subsample_columns]) ; row.names(cts) <- x$Geneid # ceiling is to round up
  # defining metadata
  coldata <- data.frame(
    condition=subsamples$condition,
    type=rep(args[5], length(subsamples$condition))
  )
  row.names(coldata) <- paste(subsamples$condition, subsamples$sample, sep="_")
  if (!all(rownames(coldata) == colnames(cts))) {
    print("the names in count_file and the sample table do not match")
  } # IF FALSE YOU ARE IN DEEP SHIT MY MAN, NOT GONNA WORK
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
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
upTEGs <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
upTEs <- res[row.names(res) %in% TEs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
upPCGs <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.05,]
downTEGs <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
downTEs <- res[row.names(res) %in% TEs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
downPCGs <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.05,]
DEGs <- list(upTEGs, upTEs, upPCGs, downTEGs, downTEs, downPCGs)
names(DEGs) <- c("upTEGs", "upTEs", "upPCGs", "downTEGs", "downTEs", "downPCGs")

# fill DEG tables with TPM values for all samples
TPM_to_df <- function(df, sample_nb) {
  df <- merge(as.data.frame(df), cts_summary_TPM, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}
DEGs <- lapply(DEGs, FUN=TPM_to_df, sample_nb=length(samples$condition))
# creates dataframe with average values among replicates
rep_mean <- function(x) { # calculate the average for every x columns (x is the number of replicates) and returns a matrix that we convert to a dataframe
  y <- as.data.frame(apply(array(as.matrix(x[sample_columns]), c(nrow(x), nb_replicate, length(unique(samples$condition)))), nb_replicate, rowMeans, na.rm = TRUE))
  names(y) <- unique(gsub(x=names(cts_summary_TPM)[sample_columns[1]:ncol(cts_summary_TPM)], pattern="_R\\d", replacement=""))
  return(y)
}
cts_summary_TPM_mean <- cbind(cts_summary_TPM[which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], rep_mean(cts_summary_TPM))

# create DEG directory
dirname <- paste("DEGs_", args[3], "_vs_", args[2], "/", sep="")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)

# plot heatmaps
suppressPackageStartupMessages(library("RColorBrewer")) ; suppressPackageStartupMessages(library("pheatmap"))
DEG_heatmap <- function(x, y, z) {
  if (is.data.frame(x)==T) {
    if (nrow(x) > 3) {
      sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x$Geneid)+1))
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(dirname, "heatmap_", y, "_", args[3], "_vs_", args[2], ".pdf"), main=paste0(args[3], "_", y, "\nn=", nrow(x),"\nlog2(tpm+1)"))
    }
  }
}
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_TPM_mean[,which(names(cts_summary_TPM_mean) %in% samples$condition)]) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_TPM[,sample_columns]) )

# export DEG tables
mapply(FUN=write.table, x=DEGs, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)

#### DEGs in antisense, filtering out DEGs that are already misregulated in the same direction in sense
cts_summary_AS <- read.delim("counts_AS_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="")
res_AS <- DESeq2_function(cts_summary_AS)
dirname <- paste("DEGs_", args[3], "_vs_", args[2], "_AS/", sep="")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
# DEGs with AS ####
upTEGs_AS <- res_AS[row.names(res_AS) %in% TEGs$GeneId & res_AS$log2FoldChange >=1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(upTEGs)),]
upTEs_AS <- res_AS[row.names(res_AS) %in% TEs$GeneId & res_AS$log2FoldChange >=1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(upTEs)),]
upPCGs_AS <- res_AS[row.names(res_AS) %in% PCGs$GeneId & res_AS$log2FoldChange >=1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(upPCGs)),]
downTEGs_AS <- res_AS[row.names(res_AS) %in% TEGs$GeneId & res_AS$log2FoldChange <=-1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(downTEGs)),]
downTEs_AS <- res_AS[row.names(res_AS) %in% TEs$GeneId & res_AS$log2FoldChange <=-1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(downTEs)),]
downPCGs_AS <- res_AS[row.names(res_AS) %in% PCGs$GeneId & res_AS$log2FoldChange <=-1 & !is.na(res_AS$padj) & res_AS$padj < 0.05 & !(row.names(res_AS) %in% row.names(downPCGs)),]
DEGs_AS <- list(upTEGs_AS, upTEs_AS, upPCGs_AS, downTEGs_AS, downTEs_AS, downPCGs_AS)
names(DEGs_AS) <- c("upTEGs_AS", "upTEs_AS", "upPCGs_AS", "downTEGs_AS", "downTEs_AS", "downPCGs_AS")
DEGs_AS <- lapply(DEGs_AS, FUN=TPM_to_df, sample_nb=length(samples$condition))
# plot heatmaps
mapply(FUN = DEG_heatmap, x=DEGs_AS, y=paste0(names(DEGs_AS), "_mean"), MoreArgs = list(z=cts_summary_TPM_mean[,which(names(cts_summary_TPM_mean) %in% samples$condition)]) )
mapply(FUN = DEG_heatmap, x=DEGs_AS, y=names(DEGs_AS), MoreArgs = list(z=cts_summary_TPM[,sample_columns]) )

# export DEG tables
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
mapply(FUN=write.table, x=DEGs_AS, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs_AS), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs_AS, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_AS_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)

unlink("Rplots.pdf") # delete that plot that is always created for some reason