#!/usr/bin/env Rscript

# ============================================================================
# DESeq2_pairwise_contrasts.R
#
# Companion to DESeq2_pipeline.R: extracts EXTRA pairwise contrasts (e.g.
# combined double-mutant vs single-mutant, rather than vs the dataset's
# WT/reference) from an ALREADY-FITTED Master dds object produced by a prior
# DESeq2_pipeline.R run, without re-estimating dispersion.
#
# Rationale: DESeq2's per-gene dispersion fit does not depend on which
# condition is chosen as the model's reference level (relevel() just
# reparametrizes the same design matrix), so results() with an explicit
# contrast= is already reference-independent and could in principle be read
# straight off the existing dds. The one exception is lfcShrink(type="apeglm"),
# which needs a named model COEFFICIENT rather than an arbitrary contrast. To
# get that coefficient for an arbitrary (target, baseline) pair without a full
# DESeq() re-run, this script relevels a COPY of the dds to the desired
# baseline and re-fits only the Wald test (nbinomWaldTest()), which reuses the
# dispersions already stored in the object instead of recalculating them --
# seconds, not minutes -- and keeps the same apeglm-shrunk log2FC methodology
# the main pipeline uses for every other comparison.
#
# Every other piece of output (DEG filtering/thresholds, TPM/RPM attachment,
# heatmaps, table layout) mirrors DESeq2_pipeline.R's per-comparison loop
# exactly, so DEGs_<combined>_vs_<single>[_AS]/ folders produced here look and
# behave identically to the DEGs_<cond>_vs_<ref>[_AS]/ folders from the
# original run -- just extra ones, addressing a different baseline.
#
# args[1] = path to the mapping-pipeline folder (read-only input; where
#           counts_summary.tsv / counts_AS_summary.tsv live) -- same as
#           DESeq2_pipeline.R's args[1] for this dataset.
# args[2] = path to the sample table used for the ORIGINAL master run (same
#           file, comma-separated, columns 2 & 3 "condition"/"sample").
# args[3] = "sense" or "antisense" -- must match the strand of the dds[_AS].rds
#           being loaded.
# args[4] = "quantseq" (RPM) or "" (TPM) -- must match the normalization used
#           in the original run (so the TPM/RPM columns attached to DEG tables
#           here are computed the same way).
# args[5] = dds_dir: READ-ONLY, the output directory of a completed
#           DESeq2_pipeline.R run for this dataset (holds dds[_AS].rds).
#           Never written to.
# args[6] = comma-separated list of "combined"/double-mutant conditions
#           (e.g. h2a.w KO + a silencing-mutant), no spaces.
# args[7] = comma-separated list of single-mutant conditions to use as the
#           baseline for each contrast, SAME LENGTH AND ORDER as args[6] --
#           args[6][i] is contrasted against args[7][i].
# args[8] = optional output directory for the new
#           DEGs_<combined>_vs_<single>[_AS]/ folders and
#           DEGs_manifest_custom[_AS].tsv. Defaults to args[5] (old behavior:
#           write alongside the original run) when omitted -- same
#           input_dir-vs-output_dir convention as DESeq2_pipeline.R's
#           args[1]/args[7], so a fresh dated run doesn't require copying
#           dds[_AS].rds anywhere first.
#
# NO HYPHEN in condition/sample names, only "." and "_" are safe (same rule as
# DESeq2_pipeline.R).
#
# Example invocation (w2_c23_s456_a56_m1_10d_in_vitro_qseq dataset: each
# silencing single mutant used as baseline for its h2a.w-combined counterpart,
# written into a fresh dated output folder alongside the original run):
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_c23_s456_a56_m1_10d_in_vitro_qseq",
#           "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_c23_s456_a56_m1_qseq.samples",
#           "sense", "quantseq",
#           "/groups/berger/user/pierre.bourguet/projects/2020_h2aw_suvh456/04_results/genomics/pipeline_runs/rnaseq/w2_c23_s456_a56_m1_10d_in_vitro_qseq/2026-08-26_DESeq2",
#           "w_cmt23,w_suvh456,w_atxr56,w_met1",
#           "cmt23,suvh456,atxr56,met1",
#           "/groups/berger/user/pierre.bourguet/projects/2020_h2aw_suvh456/04_results/genomics/pipeline_runs/rnaseq/w2_c23_s456_a56_m1_10d_in_vitro_qseq/2026-08-27_DESeq2_single_vs_combined")
# ============================================================================

# ---- Libraries ----
.libPaths(c("/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/Rlib", .libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(apeglm))

# ---- User-customizable thresholds (kept in sync BY HAND with DESeq2_pipeline.R) ----
lfc_threshold  <- 1     # |log2FoldChange| cutoff (on the un-shrunken results() estimate)
padj_threshold <- 0.05  # adjusted p-value (BH) cutoff

# ---- Command-line arguments ----
args <- commandArgs(TRUE)

input_dir <- ifelse(endsWith(args[1], "/"), args[1], paste0(args[1], "/"))
sample_table_path <- args[2]
strand_arg <- args[3]
norm_arg   <- args[4]
dds_dir <- ifelse(endsWith(args[5], "/"), args[5], paste0(args[5], "/")) # read-only
combined_conditions <- strsplit(args[6], ",")[[1]]
single_conditions   <- strsplit(args[7], ",")[[1]]
output_dir <- if (length(args) >= 8 && nzchar(args[8])) args[8] else args[5]
output_dir <- ifelse(endsWith(output_dir, "/"), output_dir, paste0(output_dir, "/"))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

if (length(combined_conditions) != length(single_conditions)) {
  stop("args[6] (combined conditions) and args[7] (single-mutant baselines) must be the same length and in matching order")
}
if (length(combined_conditions) == 0) {
  stop("no condition pairs supplied (args[6]/args[7] empty)")
}

# ============================================================================
# Annotation import (Araport11 SAF: PCGs, TEs, TEGs) -- identical to DESeq2_pipeline.R
# ============================================================================
TEGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="transposable_element_gene")
TEs  <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="transposable_element")
PCGs <- subset(read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE), Type=="gene")

# ============================================================================
# Import counts, apply sense/antisense selection (args[3]) -- identical logic
# to DESeq2_pipeline.R's args[5], just re-numbered
# ============================================================================
if (strand_arg == "sense") {
  cts_summary <- read.delim(paste0(input_dir, "counts_summary.tsv"), header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
  strand <- ""
} else if (strand_arg == "antisense") {
  cts_summary <- read.delim(paste0(input_dir, "counts_AS_summary.tsv"), header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
  strand <- "_AS"
} else {
  stop("args[3] should be 'sense' or 'antisense'")
}
setwd(output_dir)
row.names(cts_summary) <- cts_summary$Geneid
cts_summary <- cts_summary[!(substr(cts_summary$Chr,4,4) == "C" | substr(cts_summary$Chr,4,4) == "M"),] # removing chloroplastic and mitochondrial genes

samples <- read.delim(sample_table_path, header=T, sep=',', quote="", dec=".", comment.char="", check.names=FALSE)[,c(2,3)]
samples <- samples[order(samples[,1], samples[,2]),] # reorder samples by alphabetical names
sample_columns <- which(names(cts_summary) %in% paste(samples$condition, samples$sample, sep="_")) # columns holding sample counts

# Drop any genotype/replicate column present in counts_summary.tsv but not listed in the sample
# table (e.g. when input_dir is shared with a larger dataset and this run's sample table is a
# genuine subset of its genotypes) -- annotate_df()/average_replicates() below assume columns 1:8
# are metadata and 9:ncol are exactly this run's samples; leaving extra columns in would leak raw,
# unnormalized counts for the dropped genotypes into the exported tables instead of being excluded.
cts_summary <- cts_summary[, c(1:8, sample_columns)]
sample_columns <- 9:ncol(cts_summary)

missing_conditions <- setdiff(c(combined_conditions, single_conditions), unique(samples$condition))
if (length(missing_conditions) > 0) {
  stop("condition(s) not found in the sample table: ", paste(missing_conditions, collapse=", "))
}

# ============================================================================
# TPM / RPM normalization (cts_summary_norm) -- identical to DESeq2_pipeline.R,
# needed to re-derive the same per-sample + mean-of-replicates columns
# attached to DEG tables/heatmaps below (these per-sample values are not
# persisted to disk by the original run, only their replicate-mean is).
# ============================================================================
exon_size <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_gtf_exonic_genes_sizes.tsv", head=F, sep="\t", quote="", comment.char="", col.names=c("Geneid", "Length"))
exon_size_tmp <- merge(cts_summary, exon_size, by="Geneid", sort=F)
if (sum(as.character(cts_summary$Geneid[cts_summary$Geneid %in% exon_size_tmp$Geneid]) != as.character(exon_size_tmp$Geneid))>0) {
  print("re-ordering of dataframes went wrong")
}
cts_summary_norm <- cts_summary
cts_summary_norm$Length[cts_summary_norm$Geneid %in% exon_size_tmp$Geneid] <- exon_size_tmp$Length.y

RPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; }
RPK <- function(x) { ; return( x / (cts_summary_norm$Length / 1000) ) ; }
TPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; }
if (nzchar(norm_arg)) {
  if (norm_arg == "quantseq") {
    cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = RPM)
    normalization <- "RPM"
  } else {
    stop("args[4] should be 'quantseq' or left empty")
  }
} else {
  cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = RPK)
  cts_summary_norm[,sample_columns] <- apply(X=cts_summary_norm[,sample_columns], MARGIN = 2, FUN = TPM)
  normalization <- "TPM"
}

average_replicates <- function(df) {
  df <- df[order(df$Chr, df$Start),]
  df_mean <- cts_summary_norm[cts_summary_norm$Geneid %in% row.names(df), which(!names(cts_summary_norm) %in% names(cts_summary_norm)[sample_columns])]
  df_mean <- df_mean[order(df_mean$Chr, df_mean$Start),]
  for (i in unique(samples$condition)) {
    x <- df[, which(gsub(pattern=paste0("_", substr(samples$sample[1], 1, nchar( as.character(samples$sample[1])) -1), "\\d"), replacement="" , x = names(df)) %in% i)]
    if (is.data.frame(x)) {
      df_mean <- cbind(df_mean, rowMeans(x))
    } else {
      df_mean <- cbind(df_mean, x)
    }
    names(df_mean)[ncol(df_mean)] <- i
  }
  return(df_mean)
}
cts_summary_norm_mean <- average_replicates(cts_summary_norm)

annotate_df <- function(df) {
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}
mean_TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm_mean) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1]
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}
TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1]
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}

# ============================================================================
# Load the already-fitted Master dds (NOT re-fit: DESeq() -- and therefore
# dispersion estimation -- already ran once in the original DESeq2_pipeline.R
# job that produced this file).
# ============================================================================
dds_path <- paste0(dds_dir, "dds", strand, ".rds")
if (!file.exists(dds_path)) stop("could not find ", dds_path, " -- run DESeq2_pipeline.R for this dataset/strand first")
dds <- readRDS(dds_path)

missing_levels <- setdiff(c(combined_conditions, single_conditions), levels(dds$condition))
if (length(missing_levels) > 0) {
  stop("condition(s) not found among dds$condition levels (check spelling/strand): ", paste(missing_levels, collapse=", "))
}

# ============================================================================
# Heatmap helper -- identical to DEG_heatmap_df in DESeq2_pipeline.R
# ============================================================================
DEG_heatmap_df <- function(x, y, z, dirname) {
  if (is.data.frame(x)==T) {
    if (nrow(x) > 3) {
      sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x$Geneid)+1))
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(dirname, "heatmap_", y, ".pdf"), main=paste0(y, "\nn=", nrow(x),"\nlog2(", normalization, "+1)"), cluster_cols = F)
    }
  }
}

# ============================================================================
# Pairwise loop: relevel a COPY of dds to each single-mutant baseline, refit
# ONLY the Wald test (dispersions reused as-is), extract + shrink (apeglm) +
# filter + export -- same shape as DESeq2_pipeline.R's per-comparison loop.
# ============================================================================
manifest_file <- paste0("DEGs_manifest_custom", strand, ".tsv")
if (file.exists(manifest_file)) file.remove(manifest_file)

for (i in seq_along(combined_conditions)) {
  combined <- combined_conditions[i]
  single   <- single_conditions[i]
  message("Processing custom pairwise comparison: ", combined, " vs ", single)

  dds_i <- dds # cheap copy-on-write; leaves the loaded master `dds` untouched for the next iteration
  dds_i$condition <- relevel(dds_i$condition, ref = single)
  dds_i <- nbinomWaldTest(dds_i) # reuses dispersions already stored in dds_i -- does NOT recompute them

  res <- results(dds_i, contrast=c("condition", combined, single))

  coef_name <- paste0("condition_", combined, "_vs_", single)
  if (!(coef_name %in% resultsNames(dds_i))) {
    stop("expected coefficient '", coef_name, "' not found in resultsNames(dds_i) -- got: ", paste(resultsNames(dds_i), collapse=", "))
  }
  res_shrunk <- lfcShrink(dds_i, coef=coef_name, type="apeglm")

  upTEGs   <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange >= lfc_threshold  & !is.na(res$padj) & res$padj < padj_threshold,]
  upTEs    <- res[row.names(res) %in% TEs$GeneId  & res$log2FoldChange >= lfc_threshold  & !is.na(res$padj) & res$padj < padj_threshold,]
  upPCGs   <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange >= lfc_threshold  & !is.na(res$padj) & res$padj < padj_threshold,]
  downTEGs <- res[row.names(res) %in% TEGs$GeneId & res$log2FoldChange <= -lfc_threshold & !is.na(res$padj) & res$padj < padj_threshold,]
  downTEs  <- res[row.names(res) %in% TEs$GeneId  & res$log2FoldChange <= -lfc_threshold & !is.na(res$padj) & res$padj < padj_threshold,]
  downPCGs <- res[row.names(res) %in% PCGs$GeneId & res$log2FoldChange <= -lfc_threshold & !is.na(res$padj) & res$padj < padj_threshold,]
  DEGs <- list(upTEGs=upTEGs, upTEs=upTEs, upPCGs=upPCGs, downTEGs=downTEGs, downTEs=downTEs, downPCGs=downPCGs)

  DEGs_mean <- lapply(DEGs, FUN=mean_TPM_to_df, sample_nb=length(unique(sort(samples$condition))))
  DEGs      <- lapply(DEGs, FUN=TPM_to_df,      sample_nb=length(samples$condition))

  dirname_pair <- paste0("DEGs_", combined, "_vs_", single, strand, "/")
  ifelse(!dir.exists(dirname_pair), dir.create(dirname_pair), FALSE)

  df_shrunk <- data.frame(Geneid=row.names(res_shrunk), as.data.frame(res_shrunk), check.names=FALSE)
  write.table(df_shrunk, file=paste0(dirname_pair, combined, "_vs_", single, "_shrunken_log2FC_all_genes.tsv"), quote=F, sep="\t", row.names=F, col.names=T)

  # un-shrunken (raw results()) log2FC estimates for all genes -- same shape as the shrunken
  # export above, but from `res` (the value DEG calling itself uses) rather than `res_shrunk`.
  # Mirrors DESeq2_pipeline.R's per-comparison shrunken/un-shrunken pairing.
  df_raw <- data.frame(Geneid=row.names(res), as.data.frame(res), check.names=FALSE)
  write.table(df_raw, file=paste0(dirname_pair, combined, "_vs_", single, "_unshrunken_log2FC_all_genes.tsv"), quote=F, sep="\t", row.names=F, col.names=T)

  invisible(mapply(FUN = DEG_heatmap_df, x=DEGs_mean, y=paste0(combined, "_vs_", single, "_", names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,which(names(cts_summary_norm_mean) %in% samples$condition)], dirname=dirname_pair) ))
  invisible(mapply(FUN = DEG_heatmap_df, x=DEGs,      y=paste0(combined, "_vs_", single, "_", names(DEGs)),        MoreArgs = list(z=cts_summary_norm[,sample_columns], dirname=dirname_pair) ))

  invisible(mapply(FUN=write.table, x=DEGs,      file=paste0(dirname_pair, combined, "_vs_", single, "_", names(DEGs), ".tsv"),      MoreArgs=list(quote=F, sep="\t", row.names=F, col.names=T)))
  invisible(mapply(FUN=write.table, x=DEGs_mean, file=paste0(dirname_pair, combined, "_vs_", single, "_", names(DEGs), "_mean.tsv"), MoreArgs=list(quote=F, sep="\t", row.names=F, col.names=T)))
  write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste0(dirname_pair, combined, "_vs_", single, "_DEGs_summary.tsv"), quote=F, sep="\t", row.names=T, col.names=F)

  cat(paste0("DEGs_", combined, "_vs_", single, strand), "\n", file=manifest_file, append=TRUE, sep="")
}

unlink("Rplots.pdf")
