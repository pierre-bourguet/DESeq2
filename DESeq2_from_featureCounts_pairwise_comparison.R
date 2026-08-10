#!/usr/bin/env Rscript 

args <- commandArgs(TRUE);
# args[1] is path to main output folder
# args[2] is reference sample
# args[3] is treatment sample
# args[4] is path to sample table
# args[5] is paired-end or single-strand
# args[6] is optional and should be empty or "quantseq"

# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w1_h1_suvh4_cmt2_nrpd2a_10-d-old_seedlings/", "Col_WT", "nrpd2a_1", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w1_suvh4_nrpd2a.tab", "single-strand")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/iswi_ribodepletion_seedlings/", "WT", "dko", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/iswi_ribodepletion_seedlings.samples", "single-strand")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/F1_w2_suvh456_cmt23_10d_seedlings_qseq_merged_reciprocal_crosses/", "Col", "h2aw_2_suvh456", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/F1_w2_suvh456_cmt23_10d_seedlings_Qseq_ID_11680_merged_reciprocal_crosses.samples", "single-strand")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_hetero_kanno_tagseq/", "virtual_F1_10", "ddm1xCol_10", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/ddm1_kanno_virtual_F1.samples", "single-strand", "quantseq")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/MM2D_aphidicolin_x5_PE75/", "control", "Aph_10h", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/MM2D_aphidicolin_Anna.samples", "paired-end")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/cdca7_ab_ddm1_qseq_35d_RL_qseq/", "Col_0", "cdca7_ab", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/cdca7_ab_ddm1_qseq_35d_RL.samples", "single-strand", "quantseq")

# importing DFs
setwd(args[1])
TEGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="transposable_element_gene") # tried importing the updated SAF file with updated mitochondrial genome but this needs updated exon size for TPM normalization, and the exon size calculator tool I used failed with the updated mito GFF file
TEs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="transposable_element")
PCGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="gene")
cts_summary <- read.delim("counts_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
row.names(cts_summary) <- cts_summary$Geneid
cts_summary <- cts_summary[!(substr(cts_summary$Chr,4,4) == "C" | substr(cts_summary$Chr,4,4) == "M"),] # removing chloroplastic and mitochondrial genes
samples <- read.delim(args[4], header=T, sep=',', quote="", dec=".", comment.char="", check.names=FALSE)[,c(2,3)]
samples <- samples[order(samples[,1], samples[,2]),] # reorder samples by alphabetical names
subsamples <- subset(samples, samples$condition %in% c(args[2], args[3]))
subsample_columns <- which(names(cts_summary) %in% paste(subsamples$condition, subsamples$sample, sep="_")) # retrieve columns that contain subsample names
sample_columns <- which(names(cts_summary) %in% paste(samples$condition, samples$sample, sep="_")) # retrieve columns that contain sample names
# metadata <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""),  subset = !Chr %in% c("chrC", "chrM")) # not using this so far

# import PCG exon size to perform TPM normalization on cDNA length (for PCG only, TEs are normalized on full length)
exon_size <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_gtf_exonic_genes_sizes.tsv", head=F, sep="\t", quote="", comment.char="", col.names=c("Geneid", "Length"))
exon_size_tmp <- merge(cts_summary, exon_size, by="Geneid", sort=F) # input files should already be in the same order
if (sum(droplevels(cts_summary$Geneid[cts_summary$Geneid %in% exon_size_tmp$Geneid]) != droplevels(exon_size_tmp$Geneid))>0) {
  print("re-ordering of dataframes went wrong") # check that ordering is correct
}
cts_summary_norm <- cts_summary
cts_summary_norm$Length[cts_summary_norm$Geneid %in% exon_size_tmp$Geneid] <- exon_size_tmp$Length.y # replaces length by exon length for PCGs only

# TPM or RPM normalization
RPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; } # reads per million
RPK <- function(x) { ; return( x / (cts_summary_norm$Length / 1000) ) ; } 
TPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; } # transcripts per million
if (length(args) == 6) { # test if argument 6 has been provided
  if (args[6] == "quantseq") { # normalize by RPM
    cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = RPM) # normalizes all columns
    normalization <- "RPM"
  } else {
    stop("6th argument should be 'quantseq' or left empty")
  }
} else {
  cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = RPK) # normalizes all columns by exon size
  cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = TPM) # normalizes all columns by library size
  normalization <- "TPM"
} 

# creates dataframe with average values among all replicates.
cts_summary_norm_mean <- cts_summary_norm[which(!names(cts_summary) %in% names(cts_summary)[sample_columns])]
for (i in unique(samples$condition)) {
  x <- cts_summary_norm[,which(gsub(pattern="_R\\d", replacement="" , x = names(cts_summary_norm)) %in% i)]
  if (is.data.frame(x)) { # this is TRUE for most cases, when there is more than 1 replicate
    cts_summary_norm_mean <- cbind(cts_summary_norm_mean, rowMeans(x))
  } else { # this is to handle unreplicated data
    cts_summary_norm_mean <- cbind(cts_summary_norm_mean, x)
  }
  names(cts_summary_norm_mean)[ncol(cts_summary_norm_mean)] <- i
}

# DESeq2 function
suppressPackageStartupMessages(library("DESeq2"))
DESeq2_function <- function(x) { # x should be cts_summary.tsv file (summary of raw counts)
  cts <- ceiling(x[,subsample_columns]) ; row.names(cts) <- x$Geneid # ceiling is to round up
  cts <- cts[,order(names(cts))] # reorder cts by alphabetical names
  # defining metadata
  coldata <- data.frame(
    condition=subsamples$condition,
    type=rep(args[5], length(subsamples$condition))
  )
  row.names(coldata) <- paste(subsamples$condition, subsamples$sample, sep="_")
  coldata <- coldata[order(rownames(coldata)),] # reorder coldata alphabetically to match the cts column ordering above
  if (!all(rownames(coldata) == colnames(cts))) {
    stop("the names in count_file and the sample table do not match")
  }
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

# fill DEG tables with average TPM values for all samples
mean_TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm_mean) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}
DEGs_mean <- lapply(DEGs, FUN=mean_TPM_to_df, sample_nb=length(unique(sort(samples$condition))))

# fill DEG tables with TPM values for all samples
TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}
DEGs <- lapply(DEGs, FUN=TPM_to_df, sample_nb=length(samples$condition))

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
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(dirname, "heatmap_", args[3], "_vs_", args[2], "_", y, ".pdf"), main=paste0(args[3], "_", y, "\nn=", nrow(x),"\nlog2(", normalization, "+1)"), cluster_cols = F)
    }
  }
}
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,which(names(cts_summary_norm_mean) %in% samples$condition)]) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_norm[,sample_columns]) )

# export DEG tables
mapply(FUN=write.table, x=DEGs, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T)) # export tables with normalized counts
mapply(FUN=write.table, x=DEGs_mean, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs), "_mean.tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T)) # export tables with averaged normalized counts
write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F) # export summary table with number of DEGs

#### DEGs in antisense, filtering out DEGs that are already misregulated in the same direction in sense
cts_summary_AS <- read.delim("counts_AS_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
row.names(cts_summary_AS) <- cts_summary_AS$Geneid
cts_summary_AS <- cts_summary_AS[!(substr(cts_summary_AS$Chr,4,4) == "C" | substr(cts_summary_AS$Chr,4,4) == "M"),] # removing chloroplastic and mitochondrial genes, consistent with the sense counts
stopifnot("sense and antisense count files must have identical column layouts" = identical(names(cts_summary), names(cts_summary_AS))) # subsample_columns / sample_columns indices are reused for the AS table
res_AS <- DESeq2_function(cts_summary_AS)

# TPM / RPM normalization of the antisense counts (mirrors the sense normalization above)
cts_summary_AS_norm <- cts_summary_AS
cts_summary_AS_norm$Length <- cts_summary_norm$Length[match(cts_summary_AS_norm$Geneid, cts_summary_norm$Geneid)] # reuse exon-corrected lengths
if (normalization == "RPM") {
  cts_summary_AS_norm[,sample_columns] <- apply(X=cts_summary_AS_norm[,sample_columns], MARGIN = 2, FUN = RPM)
} else {
  RPK_AS <- function(x) { ; return( x / (cts_summary_AS_norm$Length / 1000) ) ; }
  cts_summary_AS_norm[,sample_columns] <- apply(X=cts_summary_AS_norm[,sample_columns], MARGIN = 2, FUN = RPK_AS)
  cts_summary_AS_norm[,sample_columns] <- apply(X=cts_summary_AS_norm[,sample_columns], MARGIN = 2, FUN = TPM)
}
# average AS values over replicates
cts_summary_AS_norm_mean <- cts_summary_AS_norm[which(!names(cts_summary_AS) %in% names(cts_summary_AS)[sample_columns])]
for (i in unique(samples$condition)) {
  x <- cts_summary_AS_norm[,which(gsub(pattern="_R\\d", replacement="" , x = names(cts_summary_AS_norm)) %in% i)]
  if (is.data.frame(x)) { # this is TRUE for most cases, when there is more than 1 replicate
    cts_summary_AS_norm_mean <- cbind(cts_summary_AS_norm_mean, rowMeans(x))
  } else { # this is to handle unreplicated data
    cts_summary_AS_norm_mean <- cbind(cts_summary_AS_norm_mean, x)
  }
  names(cts_summary_AS_norm_mean)[ncol(cts_summary_AS_norm_mean)] <- i
}
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
DEGs_AS_mean <- lapply(DEGs_AS, FUN=mean_TPM_to_df, sample_nb=length(unique(sort(samples$condition))), norm_df=cts_summary_AS_norm_mean) # must be computed BEFORE DEGs_AS is overwritten by TPM_to_df (row.names merge would break otherwise)
DEGs_AS <- lapply(DEGs_AS, FUN=TPM_to_df, sample_nb=length(samples$condition), norm_df=cts_summary_AS_norm)
# plot heatmaps (using antisense-normalized values)
mapply(FUN = DEG_heatmap, x=DEGs_AS, y=paste0(names(DEGs_AS), "_mean"), MoreArgs = list(z=cts_summary_AS_norm_mean[,which(names(cts_summary_AS_norm_mean) %in% samples$condition)]) )
mapply(FUN = DEG_heatmap, x=DEGs_AS, y=names(DEGs_AS), MoreArgs = list(z=cts_summary_AS_norm[,sample_columns]) )

# export DEG tables
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
mapply(FUN=write.table, x=DEGs_AS, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs_AS), ".tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
mapply(FUN=write.table, x=DEGs_AS_mean, file=paste(dirname, args[3], "_vs_", args[2], "_", names(DEGs_AS), "_mean.tsv", sep=""), MoreArgs=list(quote = F, sep="\t", row.names=F, col.names=T))
write.table(x=t(as.data.frame(lapply(DEGs_AS, FUN=nrow))), file=paste(dirname, args[3], "_vs_", args[2], "_AS_DEGs_summary", ".tsv", sep=""), quote = F, sep="\t", row.names=T, col.names=F)

unlink("Rplots.pdf") # delete that plot that is always created for some reason