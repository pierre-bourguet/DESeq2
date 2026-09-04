#!/usr/bin/env Rscript

args <- commandArgs(TRUE);
# args[1] is path to main output folder
# args[2] is reference sample
# args[3] is path to sample table
# args[4] is paired-end or single-strand
# args[5] is sense or antisense
# args[6] is optional and should be empty or "quantseq"

# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3_polyA/", "WT", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_polyA.samples", "paired-end", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_met1_suvh456_noWT_polyA_4-5w_leaves/", "Col0", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/ddm1_met1_suvh456_withWT_polyA.samples", "single-strand", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w1_h1_suvh4_cmt2_nrpd2a/", "WT_Col", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w1_suvh4_nrpd2a.tab", "single-strand", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/pol2a_12_polyA_13d_seedlings/", "L5", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/pol2a_12_polyA_13d_seedlings.samples", "single-strand", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/2020_Rougee_ddm1_clf/", "Col", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/2020_Rougee_ddm1_clf.samples", "single-strand", "antisense", "quantseq")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w_heat_stress_qseq/", "WT_23", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w_heat_stress_qseq.samples", "single-strand", "sense", "quantseq")
# args <- c("/groups/berger/lab/pierre.bourguet/W6_complementation/w_cmt23_screen/", "WT", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w_cmt23_W6_complementation_screen.samples", "single-strand", "sense", "quantseq")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/tagseq_01_all_reads_new_barcodes/15978/deseq2/", "Col_0", "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/tagseq01_15978_try3.csv", "single-strand", "sense", "quantseq")
# args <- c("/groups/berger/user/pierre.bourguet/projects/03_regulatory_TEs/results/RNAseq/Zhao_2019_Dev_Cell", "EC", "/groups/berger/user/pierre.bourguet/projects/03_regulatory_TEs/data/Zhao_2019/sample_sheet.csv", "paired-end", "sense")

# metadata <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""),  subset = !Chr %in% c("chrC", "chrM"))  # not using this so far
TEGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element_gene")
TEs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="transposable_element")
PCGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char=""), Type=="gene")

# importing DFs
setwd(args[1])
if (args[5]=="sense") {
  cts_summary <- read.delim("counts_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="", check.names = FALSE)
  strand <- ""
} else if (args[5]=="antisense") {
  cts_summary <- read.delim("counts_AS_summary.tsv", header=T, sep='\t', quote="", dec=".", comment.char="", check.names = FALSE)
  strand <- "_AS"
} else { print("5th argument should be sense or antisense") }
row.names(cts_summary) <- cts_summary$Geneid
cts_summary <- cts_summary[!(substr(cts_summary$Chr,4,4) == "C" | substr(cts_summary$Chr,4,4) == "M"),] # removing chloroplastic and mitochondrial genes
samples <- read.delim(args[3], header=T, sep=',', quote="", dec=".", comment.char="", check.names = FALSE)[,c(2,3)]
samples <- samples[order(samples[,1], samples[,2]),] # reorder samples by alphabetical names
sample_columns <- which(names(cts_summary) %in% paste(samples$condition, samples$sample, sep="_")) # retrieve columns that contain sample names
conditions <- unique(samples$condition)[unique(samples$condition) != args[2]]

# import PCG exon size to perform TPM normalization on cDNA length (for PCG only, TEs are normalized on full length)
exon_size <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_gtf_exonic_genes_sizes.tsv", head=F, sep="\t", quote="", comment.char="", col.names=c("Geneid", "Length"))
exon_size_tmp <- merge(cts_summary, exon_size, by="Geneid", sort=F) # input files should already be in the same order
if (sum(droplevels(cts_summary$Geneid[cts_summary$Geneid %in% exon_size_tmp$Geneid]) != droplevels(exon_size_tmp$Geneid))>0) {
  print("re-ordering of dataframes went wrong") # check that ordering is correct
}
cts_summary_norm <- cts_summary
cts_summary_norm$Length[cts_summary_norm$Geneid %in% exon_size_tmp$Geneid] <- exon_size_tmp$Length.y # replaces length by exon length

# TPM or RPM normalization
RPM <- function(x) { ; return( x / (sum(x) / 1000000) ) ; }
RPK <- function(x) { ; return( x / (cts_summary_norm$Length / 1000) ) ; }
TPM <- function(x) { ; return( x / (sum(x) / 1000000) ) ; }
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

# function that returns dataframe with average values among all replicates.
average_replicates <- function(df) {
  df <- df[order(df$Chr, df$Start),] # reorder to be sure cbind is done properly (see below)
  df_mean <- cts_summary_norm[cts_summary_norm$Geneid %in% row.names(df), which(!names(cts_summary_norm) %in% names(cts_summary_norm)[sample_columns])]
  df_mean <- df_mean[order(df_mean$Chr, df_mean$Start),] # reorder
  for (i in unique(samples$condition)) {
    x <- df[, which(gsub(pattern=paste0("_", substr(samples$sample[1], 1, nchar( as.character(samples$sample[1])) -1), "\\d"), replacement="" , x = names(df)) %in% i)]
    if (is.data.frame(x)) { # this is TRUE for most cases, when there is more than 1 replicate
      df_mean <- cbind(df_mean, rowMeans(x))
    } else { # this is to handle unreplicated data
      df_mean <- cbind(df_mean, x)
    }
    names(df_mean)[ncol(df_mean)] <- i
  }
  return(df_mean)
}
cts_summary_norm_mean <- average_replicates(cts_summary_norm)

# DESeq2 function ####
suppressPackageStartupMessages(library("DESeq2"))
DESeq2_function <- function(x) { # x should be cts_summary.tsv file (summary of raw counts)
  cts <- ceiling(x[,sample_columns]) ; row.names(cts) <- x$Geneid # ceiling is to round up
  # defining metadata
  coldata <- data.frame(
    condition=substr(names(cts_summary)[sample_columns], 1, nchar(names(cts_summary)[sample_columns])-3),
    type=rep(args[4], nrow(samples))
  )
  row.names(coldata) <- names(cts_summary)[sample_columns]
  if (!all(rownames(coldata) == colnames(cts))) {
    stop("the names in count_file and the sample table do not match")
  }
  dds <<- DESeqDataSetFromMatrix(countData = cts,
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
  return(DESeq(dds))
}
dds <- DESeq2_function(cts_summary)

# create output directory ####
output_dir <- paste0("DEGs_batch", strand, "/")
ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
# plotting PCA ####
pdf(paste0(output_dir, "PCA.pdf"))
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
graphics.off()
# plotting PCA on transposable elements ####
# pdf(paste0(output_dir, "PCA_TE.pdf"))
# vsd_TE <- vst(DESeq2_function(cts_summary[cts_summary$Type=="transposable_element",]), blind=FALSE, nsub=100) # you might have to lower nsub to lower than default (1000) if there are too few reads at TEs
# plotPCA(vsd_TE, intgroup=c("condition"))
# graphics.off()
# plotting PCA on transposable element genes
# pdf(paste0(output_dir, "PCA_TEG.pdf"))
# vsd_TEG <- vst(DESeq2_function(cts_summary[cts_summary$Type=="transposable_element_gene",]), blind=FALSE, nsub=100) # you might have to lower nsub to lower than default (1000) if there are too few reads at TEs
# plotPCA(vsd_TEG, intgroup=c("condition"))
# graphics.off()

# compare all mutants with control ####
f <- function(aa, bb) { # this creates a res_mutant dataframe comparing mutant vs control
  eval(substitute( a <- results(dds, contrast=c("condition",as.character(b),args[2]))
                  , list(a = aa, b = bb))) # this second argument provides an environment for substitute, defining variables used in previous line
}
all_res <- Map(f, paste0("res_", conditions), as.list(conditions)) # list of all mutant vs WT comparisons

# creates list of all up & down TEs, TEGs, PCGs (across all mutant conditions) ####
up <- function(x, y) {
  return( row.names(x[row.names(x) %in% eval(parse(text=paste0(y, "$GeneId"))) & x$log2FoldChange >=1 & !is.na(x$padj) & x$padj < 0.05,]) )
}
down <- function(x, y) {
  return( row.names(x[row.names(x) %in% eval(parse(text=paste0(y, "$GeneId"))) & x$log2FoldChange <=-1 & !is.na(x$padj) & x$padj < 0.05,]) )
}
DEGs <- list(
  upTEGs = unique(unlist(lapply(FUN = up, X = all_res, y="TEGs"))),
  upTEs = unique(unlist(lapply(FUN = up, X = all_res, y="TEs"))),
  upPCGs = unique(unlist(lapply(FUN = up, X = all_res, y="PCGs"))),
  downTEGs = unique(unlist(lapply(FUN = down, X = all_res, y="TEGs"))),
  downTEs = unique(unlist(lapply(FUN = down, X = all_res, y="TEs"))),
  downPCGs = unique(unlist(lapply(FUN = down, X = all_res, y="PCGs")))
)

# plot heatmaps ####
suppressPackageStartupMessages(library("RColorBrewer")) ; suppressPackageStartupMessages(library("pheatmap"))
DEG_heatmap <- function(x, y, z, n) {
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      if (n %in% c("rlog", "vst")) {
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\n", n)
        }
      else {
        sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x)+1))
        title <- paste0(y, "\nn=", length(x),"\nlog2(", n, "+1)")
        }
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(output_dir, n, "/heatmap_", y, ".pdf"), main=title, cluster_cols = F)
    }
  }
}

mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,9:(9+length(conditions))], n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_norm[,9:(8+nrow(samples))], n = normalization) )

# write DEG tables ####
DEG_write <- function(x, y) {
  write.table(subset(cts_summary_norm_mean, subset=row.names(cts_summary_norm_mean) %in% x), file=paste0(output_dir, "batch_", y, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
}
mapply(FUN = DEG_write, x=DEGs, y=names(DEGs))
write.table(x=t(as.data.frame(lapply(DEGs, FUN=length))), file=paste0(output_dir, "DEGs_batch.tsv"), quote = F, sep="\t", row.names=T, col.names=F)

# write file with TPM / RPM values for all samples, averaged over replicates ####
write.table(x=cts_summary_norm_mean, file=paste0(args[1], normalization, strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# function that adds annotations to a df using its row.names as Geneid
annotate_df <- function(df) {
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}
# median of ratios (from DESeq2, median of ratios to geometric mean). Write tables and heatmaps ####
cts_MoR <- annotate_df( as.data.frame(counts(dds, normalized=T)) )
cts_MoR_mean <- average_replicates(cts_MoR)
write.table(x=cts_MoR, file=paste0(args[1], "ESF", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=cts_MoR_mean, file=paste0(args[1], "ESF", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_MoR[,9:(8+nrow(samples))], n = "MoR") )
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_MoR_mean[,9:(9+length(conditions))], n = "MoR") )
# rlog normalization ####
rld <- rlog(dds, blind=FALSE)
rld_df <- annotate_df( as.data.frame(assay(rld)) )
rld_df_mean <- average_replicates(rld_df)
write.table(x=rld_df, file=paste0(args[1], "rlog", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=rld_df_mean, file=paste0(args[1], "rlog", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=rld_df_mean[,9:(9+length(conditions))], n = "rlog") )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=rld_df[,9:(8+nrow(samples))], n = "rlog") )

# VST normalization (variance stabilizing transformation) ####
# vsd was already computed above (blind=FALSE) for the PCA plot, reusing it here
vsd_df <- annotate_df( as.data.frame(assay(vsd)) )
vsd_df_mean <- average_replicates(vsd_df)
write.table(x=vsd_df, file=paste0(args[1], "vst", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=vsd_df_mean, file=paste0(args[1], "vst", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=vsd_df_mean[,9:(9+length(conditions))], n = "vst") )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=vsd_df[,9:(8+nrow(samples))], n = "vst") )

unlink("Rplots.pdf") # delete that plot that is always created for some reason

# exploring differences between median of ratios normalization by DESeq and RPM / TPM ####
if (FALSE==TRUE) { # just protecting this code so it's not executed when i run the script
  test <- cts_summary[,sample_columns]
  par(pty = "s")
  plot(x=colSums(test) / colSums(test)[6], y=sizeFactors(dds), xlim=c(0.5,1.7), ylim=c(0.5,1.7), xlab=normalization)
  abline(a = 0, b=1)
  text(x=colSums(test) / colSums(test)[6], sizeFactors(dds), labels=names(sizeFactors(dds)), cex= 0.5, pos=1)
}