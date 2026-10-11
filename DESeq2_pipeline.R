#!/usr/bin/env Rscript

# ============================================================================
# DESeq2_pipeline.R
#
# Unified DESeq2 pipeline (merges the old *_batch_analysis.R and
# *_pairwise_comparison.R scripts).
#
# Rationale for the merge: instead of building a small, separate dds object
# per Mutant-vs-WT pair (subset to just those 2 conditions), we now fit ONE
# "Master" dds across ALL samples in the experiment (design = ~ condition).
# This lets DESeq2 estimate dispersion from the full set of replicates,
# which matters most for low-count features (e.g. transposable elements).
# All standard pairwise contrasts (Mutant vs Reference) are then extracted
# from that single Master model in a loop.
#
# args[1] = path to the mapping-pipeline folder (read-only input; where
#           counts_summary.tsv / counts_AS_summary.tsv live)
# args[2] = reference / WT condition (must match a value in the sample table)
# args[3] = path to sample table (comma-separated; columns 2 & 3 must be
#           "condition" and "sample"; replicates must be suffixed R1, R2, R3...)
# args[4] = "single-strand" or "paired-end"
# args[5] = "sense" or "antisense" -> selects counts_summary.tsv vs
#           counts_AS_summary.tsv, and suffixes output with "_AS" for antisense
# args[6] = optional, "quantseq" (RPM normalization) or leave empty (TPM)
# args[7] = optional output directory for all DESeq2 outputs. Defaults to
#           args[1] if omitted (old behavior, mapping and DESeq2 outputs
#           mixed together). Pass a dedicated folder to keep mapping-pipeline
#           outputs and DESeq2 outputs fully separate -- this removes the old
#           need to `cp` mapping results into a fresh dated folder before
#           (re-)running DESeq2: just point args[1] at the (untouched,
#           read-only) mapping folder and args[7] at wherever you want this
#           run's results to go. args[1] itself is never written to.
#
# NO HYPHEN in condition/sample names, only "." and "_" are safe.
#
# The script is safe to run twice back-to-back (once per strand): every
# output file/folder it writes is suffixed by `strand` ("" for sense,
# "_AS" for antisense), so a sense run and an antisense run never collide,
# and re-running either one simply overwrites its own outputs.
#
# Antisense run: needs the sense run's dds.rds in the output directory (run
# sense first). Raw antisense counts, own dispersion trend and BH, size
# factors of the sense fit; antisense DEG calls are kept only if they pass the
# leakage and opposite-strand partner rules of antisense_calls.R (this
# folder). Writes antisense_leakage.tsv (leak rate per library) and adds the
# two values the rules use to the all-genes tables of each comparison. It
# also writes DEGs_both_orientations/: sense and antisense calls together,
# one entry per feature (a feature called in both orientations is kept in
# the orientation with the larger change; antisense_calls.R).
#
# The fitted Master dds object is saved as dds[_AS].rds at the root of the
# output directory (readRDS() it back to build custom contrasts without
# refitting).
#
# Manifest: at the end of the pairwise loop, the script (re)writes
# "DEGs_manifest.tsv" (or "DEGs_manifest_AS.tsv" for antisense) in the
# output directory, listing exactly the DEGs_<condition>_vs_<ref>[_AS]
# folders THIS run produced. Downstream aggregation (post_processing.sh /
# intersect_DEG_batch.sh) reads that manifest instead of globbing
# "DEGs_*" from the directory, so stale folders from a previous run (old
# thresholds, a removed condition, leftover test output, ...) can never
# silently corrupt the aggregated tables. The manifest for a given strand
# is truncated at the start of that strand's run, so a rerun always
# reflects only what actually succeeded this time.
#
# Example invocations (adapted from the old scripts; treatment is no longer
# a positional argument -- ALL non-reference conditions are processed in a
# single run):
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3_polyA/", "WT", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_polyA.samples", "paired-end", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_met1_suvh456_noWT_polyA_4-5w_leaves/", "Col0", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/ddm1_met1_suvh456_withWT_polyA.samples", "single-strand", "sense")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/2020_Rougee_ddm1_clf/", "Col", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/2020_Rougee_ddm1_clf.samples", "single-strand", "antisense", "quantseq")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_met1_suvh456_noWT_polyA_4-5w_leaves/", "Col0", "/groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/ddm1_met1_suvh456_withWT_polyA.samples", "single-strand", "sense", "", "/groups/berger/lab/pierre.bourguet/2026_h2aw_suvh456/04_results/genomics/pipeline_runs/rnaseq/ddm1_met1_suvh456/2026-08-10_DESeq2")
# ============================================================================

# ---- Libraries ----
# apeglm is not bundled in any r-bundle-bioconductor module on this cluster,
# so it's installed in a project-local library (see install_apeglm.R in this
# folder) -- prepend it to .libPaths() before loading it.
.libPaths(c("/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/Rlib", .libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(circlize)) # for colorRamp2, used by col_log2FC below
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "antisense_calls.R"))

# ---- User-customizable thresholds (used everywhere DEGs are filtered) ----
lfc_threshold  <- 1     # |log2FoldChange| cutoff (on the un-shrunken results() estimate, as in the old scripts)
padj_threshold <- 0.05  # adjusted p-value (BH) cutoff

# Diverging color scale for shrunken-log2FC heatmaps. Kept in sync by hand with
# shared/env/env_shared.R's col_log2FC (same name, same definition) rather than
# sourcing that file directly: env_shared.R does library(tidyverse), and this
# cluster's R module is missing tidyverse's purrr/stringr, so sourcing it here
# would crash the whole pipeline. Domain is [-3, 0, 3]; since real log2FC values
# rarely land exactly in that range, every heatmap that uses this rescales its
# breaks to this function's [-3, 3] domain to match its own actual data spread
# (see the log2FC branch of DEG_heatmap_ids below).
col_log2FC <- colorRamp2(c(-3, 0, 3), c("cornflowerblue", "#FFFFF0", "brown1"))

# ---- Command-line arguments ----
args <- commandArgs(TRUE)

# input_dir is read-only (mapping-pipeline output); output_dir is where this
# script writes everything. output_dir defaults to input_dir when args[7] is
# not supplied, preserving the old (mapping+DESeq2 in one folder) behavior.
input_dir  <- ifelse(endsWith(args[1], "/"), args[1], paste0(args[1], "/"))
output_dir <- if (length(args) >= 7 && nzchar(args[7])) args[7] else args[1]
output_dir <- ifelse(endsWith(output_dir, "/"), output_dir, paste0(output_dir, "/"))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# Annotation import (Araport11 SAF: PCGs, TEs, TEGs)
# ============================================================================
SAF  <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", head=T, sep="\t", quote="", comment.char="", check.names=FALSE)
TEGs <- subset(SAF, Type=="transposable_element_gene")
TEs  <- subset(SAF, Type=="transposable_element")
PCGs <- subset(SAF, Type=="gene")

# ============================================================================
# Import counts, apply sense/antisense selection (args[5])
# ============================================================================
# counts are read from input_dir; everything else is written relative to
# output_dir (set as the working directory once, below)
if (args[5] == "sense") {
  cts_summary <- read.delim(paste0(input_dir, "counts_summary.tsv"), header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
  strand <- ""
} else if (args[5] == "antisense") {
  cts_summary <- read.delim(paste0(input_dir, "counts_AS_summary.tsv"), header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
  strand <- "_AS"
} else {
  stop("5th argument should be 'sense' or 'antisense'")
}
setwd(output_dir)
row.names(cts_summary) <- cts_summary$Geneid
cts_summary <- cts_summary[!(substr(cts_summary$Chr,4,4) == "C" | substr(cts_summary$Chr,4,4) == "M"),] # removing chloroplastic and mitochondrial genes

samples <- read.delim(args[3], header=T, sep=',', quote="", dec=".", comment.char="", check.names=FALSE)[,c(2,3)]
samples <- samples[order(samples[,1], samples[,2]),] # reorder samples by alphabetical names
sample_columns <- which(names(cts_summary) %in% paste(samples$condition, samples$sample, sep="_")) # columns holding sample counts
conditions <- unique(samples$condition)[unique(samples$condition) != args[2]] # all non-reference conditions -> pairwise loop targets

# Drop any genotype/replicate column present in counts_summary.tsv but not listed in the sample
# table (e.g. when input_dir is shared with a larger dataset and this run's sample table is a
# genuine subset of its genotypes) -- every downstream annotation/averaging step (annotate_df(),
# average_replicates(), the batch-heatmap z= argument) assumes columns 1:8 are metadata and 9:ncol
# are exactly this run's samples; leaving extra columns in would leak raw, unnormalized counts for
# the dropped genotypes into vst/rlog/ESF/TPM/log2FC exports instead of being excluded.
cts_summary <- cts_summary[, c(1:8, sample_columns)]
sample_columns <- 9:ncol(cts_summary)

# ============================================================================
# TPM / RPM normalization (cts_summary_norm) -- unchanged from the old scripts
# ============================================================================
# import PCG exon size to perform TPM normalization on cDNA length (for PCG only, TEs are normalized on full length)
exon_size <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_gtf_exonic_genes_sizes.tsv", head=F, sep="\t", quote="", comment.char="", col.names=c("Geneid", "Length"))
exon_size_tmp <- merge(cts_summary, exon_size, by="Geneid", sort=F) # input files should already be in the same order
# as.character() rather than the original droplevels(): read.delim() defaulted
# to stringsAsFactors=TRUE under the R versions the old scripts were written
# for (Geneid was a factor), but on this cluster's current R (4.4.1, via
# r-bundle-bioconductor/3.19) it defaults to FALSE, so Geneid is already a
# plain character vector and droplevels() has no applicable method for it.
if (sum(as.character(cts_summary$Geneid[cts_summary$Geneid %in% exon_size_tmp$Geneid]) != as.character(exon_size_tmp$Geneid))>0) {
  print("re-ordering of dataframes went wrong") # check that ordering is correct
}
cts_summary_norm <- cts_summary
cts_summary_norm$Length[cts_summary_norm$Geneid %in% exon_size_tmp$Geneid] <- exon_size_tmp$Length.y # replaces length by exon length for PCGs only

RPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; } # reads per million
RPK <- function(x) { ; return( x / (cts_summary_norm$Length / 1000) ) ; }
TPM <- function(x) { ; return( x / (sum(x) / 1e6) ) ; } # transcripts per million
if (length(args) >= 6 && nzchar(args[6])) { # test if argument 6 has been provided and is non-empty
  # (args[6] may be an empty-string placeholder rather than truly absent,
  # e.g. when the sbatch wrapper passes args[7]/output_dir but no quantseq flag)
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

# function that returns a dataframe with average values among all replicates (generic replicate-suffix stripping, from the batch script)
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

# function that adds annotation columns back onto a df using its row.names as Geneid (from the batch script)
annotate_df <- function(df) {
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}

# functions that attach TPM/RPM-normalized data columns onto DEG tables (from the pairwise script, unchanged column arithmetic)
mean_TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm_mean) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}
TPM_to_df <- function(df, sample_nb, norm_df=cts_summary_norm) {
  df <- merge(as.data.frame(df), norm_df, by="row.names")[,-1] # add annotations to DEGs, remove first column which is row.names
  return(df[,c(8:13,7,14,1:6,15:(14+sample_nb))])
}

# ============================================================================
# Master DESeq2 model -- ALL samples together, design = ~ condition
# ============================================================================
cts <- ceiling(cts_summary[, sample_columns]) ; row.names(cts) <- cts_summary$Geneid # ceiling is to round up
coldata <- data.frame(
  condition = substr(names(cts_summary)[sample_columns], 1, nchar(names(cts_summary)[sample_columns]) - 3),
  type = rep(args[4], nrow(samples))
)
row.names(coldata) <- names(cts_summary)[sample_columns]
if (!all(rownames(coldata) == colnames(cts))) {
  stop("the names in count_file and the sample table do not match")
}
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
mcols(dds) <- DataFrame(mcols(dds))
keep <- rowSums(counts(dds)) >= 10 # pre-filtering low-count features
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = args[2]) # setting the reference/WT condition
if (strand == "_AS") { # antisense: size factors of the sense fit
  if (!file.exists("dds.rds")) stop("antisense run needs the sense fit (dds.rds) in ", output_dir, ": run the sense strand first")
  dds_sense <- readRDS("dds.rds")
  sizeFactors(dds) <- sizeFactors(dds_sense)[colnames(dds)]
}
dds <- DESeq(dds) # single dispersion estimation across the whole experiment

if (strand == "_AS") {
  # strand leakage per library, from the raw counts before the pre-filter
  S_tab <- read.delim(paste0(input_dir, "counts_summary.tsv"), header=T, sep='\t', quote="", dec=".", comment.char="", check.names=FALSE)
  S_raw <- as.matrix(ceiling(S_tab[, colnames(cts)])) ; row.names(S_raw) <- S_tab$Geneid
  A_raw <- as.matrix(cts)
  leak <- leak_rates(S_raw, A_raw, sizeFactors(dds), PCGs$GeneId)
  write.table(data.frame(sample=colnames(dds), condition=dds$condition, size_factor=sizeFactors(dds), leak_rate=leak[colnames(dds)]),
              file="antisense_leakage.tsv", quote=F, sep="\t", row.names=F, col.names=T)
  # features of either fit overlapping on opposite strands
  ann <- SAF[SAF$GeneId %in% union(row.names(dds), row.names(dds_sense)), ]
  strand_pairs <- opposite_strand_pairs(data.frame(id=ann$GeneId, chr=ann$Chr, start=ann$Start, end=ann$End, strand=ann$Strand))
}

# DEG calls of one comparison; antisense calls also need to pass the rules of antisense_calls.R (`kept`)
is_call <- function(res, ann, direction, kept=rep(TRUE, nrow(res))) {
  lfc_ok <- if (direction == "up") res$log2FoldChange >= lfc_threshold else res$log2FoldChange <= -lfc_threshold
  row.names(res) %in% ann$GeneId & lfc_ok & !is.na(res$padj) & res$padj < padj_threshold & kept
}
# values of the antisense rules for one comparison (cond vs ref)
antisense_rule_values <- function(res, cond, ref) {
  res_sense <- results(dds_sense, contrast=c("condition", cond, ref))
  sense_up   <- row.names(res_sense)[is_call(res_sense, SAF, "up")]
  sense_down <- row.names(res_sense)[is_call(res_sense, SAF, "down")]
  ids <- row.names(res)
  data.frame(antisense_over_leakage = antisense_over_leakage(S_raw, A_raw, sizeFactors(dds), leak, ids, res$log2FoldChange,
                                                             colnames(dds)[dds$condition == cond], colnames(dds)[dds$condition == ref]),
             opposite_strand_sense_DEG = opposite_strand_sense_calls(ids, res$log2FoldChange, strand_pairs, sense_up, sense_down))
}

# save the fitted Master dds object itself (readRDS() it back later for custom
# contrasts -- e.g. dds <- readRDS("dds.rds") -- without re-running DESeq())
saveRDS(dds, file=paste0("dds", strand, ".rds"))

# ============================================================================
# Batch exports: PCA, VST, rlog, median-of-ratios (MoR) -- from the old batch script
# ============================================================================
batch_dir <- paste0("DEGs_batch", strand, "/")
ifelse(!dir.exists(batch_dir), dir.create(batch_dir), FALSE)

# PCA -- general (all annotation types together, DESeq2's default top-500-most-variable), plus one per annotation type (PCGs / TEGs / TEs) subset from the same vst values
vsd <- vst(dds, blind=FALSE)

pdf(paste0(batch_dir, "PCA.pdf"))
print(plotPCA(vsd, intgroup=c("condition")))
graphics.off()

for (ann_name in c("PCG", "TEG", "TE")) {
  ann_ids <- switch(ann_name, PCG=PCGs$GeneId, TEG=TEGs$GeneId, TE=TEs$GeneId)
  vsd_sub <- vsd[row.names(vsd) %in% ann_ids, ]
  if (nrow(vsd_sub) > 1) { # plotPCA needs at least 2 features to compute variance
    pdf(paste0(batch_dir, "PCA_", ann_name, ".pdf"))
    print(plotPCA(vsd_sub, intgroup=c("condition"))) # print() needed: auto-print doesn't fire for a ggplot returned from inside a for loop
    graphics.off()
  }
}

# median of ratios (DESeq2 size-factor normalized counts)
cts_MoR <- annotate_df( as.data.frame(counts(dds, normalized=T)) )
cts_MoR_mean <- average_replicates(cts_MoR)
write.table(x=cts_MoR, file=paste0("ESF", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=cts_MoR_mean, file=paste0("ESF", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# rlog
rld <- rlog(dds, blind=FALSE)
rld_df <- annotate_df( as.data.frame(assay(rld)) )
rld_df_mean <- average_replicates(rld_df)
write.table(x=rld_df, file=paste0("rlog", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=rld_df_mean, file=paste0("rlog", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# VST (vsd computed above for the PCA plot, reused here)
vsd_df <- annotate_df( as.data.frame(assay(vsd)) )
vsd_df_mean <- average_replicates(vsd_df)
write.table(x=vsd_df, file=paste0("vst", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=vsd_df_mean, file=paste0("vst", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# TPM/RPM table, averaged over replicates
write.table(x=cts_summary_norm_mean, file=paste0(normalization, strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# ============================================================================
# Heatmap helpers
#   - DEG_heatmap_ids : x is a character vector of Geneids   (used for the combined batch DEG summary below)
#   - DEG_heatmap_df  : x is a DEG dataframe with a $Geneid column (used per pairwise comparison)
# ============================================================================
DEG_heatmap_ids <- function(x, y, z, n) {
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      breaks <- NA # pheatmap default: auto breaks matching the color vector length
      if (n == "log2FC") {
        # values are already shrunken log2FC -- no rlog/vst-style scaling and no log2(x+1)
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\nshrunken log2FC")
        # col_log2FC's domain is [-3,0,3]; rescale it to the range that covers
        # the central 90% of |log2FC| here (not the full min/max) so the
        # scale isn't stretched pale by a handful of extreme outliers -- the
        # remaining 10% are winsorized (clamped) to +-lim so they still show
        # as fully saturated instead of falling outside pheatmap's breaks and
        # being drawn as NA/grey.
        lim <- suppressWarnings(as.numeric(quantile(abs(sampleDistMatrix), 0.90, na.rm=TRUE)))
        if (!is.finite(lim) || lim == 0) lim <- 1
        sampleDistMatrix[sampleDistMatrix >  lim] <-  lim
        sampleDistMatrix[sampleDistMatrix < -lim] <- -lim
        colors <- col_log2FC(seq(-3, 3, length.out=255))
        breaks <- seq(-lim, lim, length.out=256)
      } else if (n %in% c("rlog", "vst")) {
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\n", n)
        colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      } else {
        sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x)+1))
        title <- paste0(y, "\nn=", length(x),"\nlog2(", n, "+1)")
        colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      }
      rownames(sampleDistMatrix) <- NULL
      ifelse(!dir.exists(paste0(batch_dir, n)), dir.create(paste0(batch_dir, n)), FALSE)
      pheatmap(sampleDistMatrix, col=colors, breaks=breaks, filename=paste0(batch_dir, n, "/heatmap_", y, ".pdf"), main=title, cluster_cols = F)
    }
  }
}
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
# Generic pairwise extraction: loop over every non-reference condition
# Each iteration: results() + lfcShrink(normal) + DEG filtering + TPM
# attachment + heatmaps + table export ("DEGs_<condition>_vs_<ref>[_AS]/")
# The un-shrunken `res` objects are also kept (all_res) to build the
# combined cross-mutant DEG summary further below (legacy "batch" feature).
# ============================================================================
all_res <- list()
all_kept <- list() # per comparison, TRUE for features whose calls are kept (all sense features; antisense features passing the rules)
all_shrunk <- list() # shrunken-log2FC results per condition, reused below for the combined all-comparisons table and the DEGs_batch log2FC heatmap

# manifest of DEGs_<condition>_vs_<ref> folders produced by THIS run, consumed
# by post_processing.sh / intersect_DEG_batch.sh instead of them globbing
# "DEGs_*" (which would also pick up stale folders from earlier runs).
# Truncated up front so a rerun's manifest reflects only what succeeds now.
manifest_file <- paste0("DEGs_manifest", strand, ".tsv")
if (file.exists(manifest_file)) file.remove(manifest_file)

for (cond in conditions) {
  message("Processing pairwise comparison: ", cond, " vs ", args[2])

  res <- results(dds, contrast=c("condition", cond, args[2]))
  all_res[[cond]] <- res
  if (strand == "_AS") {
    rule_values <- antisense_rule_values(res, cond, args[2])
    kept <- antisense_call_kept(rule_values$antisense_over_leakage, rule_values$opposite_strand_sense_DEG)
  } else {
    kept <- rep(TRUE, nrow(res))
  }
  all_kept[[cond]] <- kept

  # shrink log2FC estimates (type="apeglm") for the whole gene set (not just DEGs). 
  coef_name <- resultsNames(dds)[startsWith(resultsNames(dds), paste0("condition_", cond, "_vs_"))]
  res_shrunk <- lfcShrink(dds, coef=coef_name, type="apeglm")
  all_shrunk[[cond]] <- res_shrunk

  # DEGs, using the customizable thresholds defined at the top of the script
  upTEGs   <- res[is_call(res, TEGs, "up", kept),]
  upTEs    <- res[is_call(res, TEs,  "up", kept),]
  upPCGs   <- res[is_call(res, PCGs, "up", kept),]
  downTEGs <- res[is_call(res, TEGs, "down", kept),]
  downTEs  <- res[is_call(res, TEs,  "down", kept),]
  downPCGs <- res[is_call(res, PCGs, "down", kept),]
  DEGs <- list(upTEGs=upTEGs, upTEs=upTEs, upPCGs=upPCGs, downTEGs=downTEGs, downTEs=downTEs, downPCGs=downPCGs)

  # attach TPM/RPM columns (mean-of-replicates version, and per-sample version)
  DEGs_mean <- lapply(DEGs, FUN=mean_TPM_to_df, sample_nb=length(unique(sort(samples$condition))))
  DEGs      <- lapply(DEGs, FUN=TPM_to_df,      sample_nb=length(samples$condition))

  # output directory for this comparison (strand-suffixed, e.g. "_AS")
  dirname_pair <- paste0("DEGs_", cond, "_vs_", args[2], strand, "/")
  ifelse(!dir.exists(dirname_pair), dir.create(dirname_pair), FALSE)

  # export shrunken log2FC estimates for all genes
  df_shrunk <- data.frame(Geneid=row.names(res_shrunk), as.data.frame(res_shrunk), check.names=FALSE)
  if (strand == "_AS") df_shrunk <- cbind(df_shrunk, rule_values)
  write.table(df_shrunk, file=paste0(dirname_pair, cond, "_vs_", args[2], "_shrunken_log2FC_all_genes.tsv"), quote=F, sep="\t", row.names=F, col.names=T)

  # heatmaps
  invisible(mapply(FUN = DEG_heatmap_df, x=DEGs_mean, y=paste0(cond, "_vs_", args[2], "_", names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,which(names(cts_summary_norm_mean) %in% samples$condition)], dirname=dirname_pair) ))
  invisible(mapply(FUN = DEG_heatmap_df, x=DEGs,      y=paste0(cond, "_vs_", args[2], "_", names(DEGs)),        MoreArgs = list(z=cts_summary_norm[,sample_columns], dirname=dirname_pair) ))

  # export DEG tables
  invisible(mapply(FUN=write.table, x=DEGs,      file=paste0(dirname_pair, cond, "_vs_", args[2], "_", names(DEGs), ".tsv"),      MoreArgs=list(quote=F, sep="\t", row.names=F, col.names=T)))
  invisible(mapply(FUN=write.table, x=DEGs_mean, file=paste0(dirname_pair, cond, "_vs_", args[2], "_", names(DEGs), "_mean.tsv"), MoreArgs=list(quote=F, sep="\t", row.names=F, col.names=T)))
  write.table(x=t(as.data.frame(lapply(DEGs, FUN=nrow))), file=paste0(dirname_pair, cond, "_vs_", args[2], "_DEGs_summary.tsv"), quote=F, sep="\t", row.names=T, col.names=F)

  # record this comparison as successfully completed (folder name, no trailing slash)
  cat(paste0("DEGs_", cond, "_vs_", args[2], strand), "\n", file=manifest_file, append=TRUE, sep="")
}

# ============================================================================
# Combined shrunken-log2FC table (all comparisons, all annotations) -- one
# column per condition-vs-reference comparison, reusing `all_shrunk` computed
# in the loop above (no re-shrinking). Exported at the root of output_dir.
# ============================================================================
log2FC_matrix <- sapply(all_shrunk, FUN = function(x) x$log2FoldChange)
row.names(log2FC_matrix) <- row.names(all_shrunk[[1]]) # same Master dds for every condition -> identical row order/IDs throughout
colnames(log2FC_matrix) <- paste0(names(all_shrunk), "_vs_", args[2])
log2FC_df <- annotate_df(as.data.frame(log2FC_matrix))
write.table(x=log2FC_df, file=paste0("shrunken_log2FC_all_comparisons", strand, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)

# ============================================================================
# Combined un-shrunken log2FC table (all comparisons, all genes) -- same
# shape as the shrunken table above, but from the raw results() estimate
# (all_res, the value DEG calling itself uses -- see lfc_threshold's comment
# at the top of this script) rather than the apeglm-shrunk one. Exported at
# the root of output_dir, transcriptome-wide (every gene that survived the
# pre-filtering step, not just DEGs).
# ============================================================================
log2FC_raw_matrix <- sapply(all_res, FUN = function(x) x$log2FoldChange)
row.names(log2FC_raw_matrix) <- row.names(all_res[[1]]) # same Master dds for every condition -> identical row order/IDs throughout
colnames(log2FC_raw_matrix) <- paste0(names(all_res), "_vs_", args[2])
log2FC_raw_df <- annotate_df(as.data.frame(log2FC_raw_matrix))
write.table(x=log2FC_raw_df, file=paste0("unshrunken_log2FC_all_comparisons", strand, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)

# ============================================================================
# Combined cross-mutant DEG summary (legacy "batch" feature): union of DEGs
# across ALL mutant conditions, reusing the `all_res` list computed above
# (no re-fitting / re-testing -- same Master model, just re-summarized).
# ============================================================================
union_calls <- function(ann, direction) {
  unique(unlist(lapply(names(all_res), function(cond) row.names(all_res[[cond]])[is_call(all_res[[cond]], ann, direction, all_kept[[cond]])])))
}
DEGs_batch <- list(
  upTEGs   = union_calls(TEGs, "up"),
  upTEs    = union_calls(TEs,  "up"),
  upPCGs   = union_calls(PCGs, "up"),
  downTEGs = union_calls(TEGs, "down"),
  downTEs  = union_calls(TEs,  "down"),
  downPCGs = union_calls(PCGs, "down")
)

# heatmaps (mean-of-replicates and per-sample, for each of TPM/RPM, MoR and rlog/vst normalizations)
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=paste0(names(DEGs_batch), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,9:(9+length(conditions))], n = normalization) ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=names(DEGs_batch),                  MoreArgs = list(z=cts_summary_norm[,9:(8+nrow(samples))],           n = normalization) ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=names(DEGs_batch),                  MoreArgs = list(z=cts_MoR[,9:(8+nrow(samples))],                    n = "MoR") ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=paste0(names(DEGs_batch), "_mean"), MoreArgs = list(z=cts_MoR_mean[,9:(9+length(conditions))],          n = "MoR") ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=paste0(names(DEGs_batch), "_mean"), MoreArgs = list(z=rld_df_mean[,9:(9+length(conditions))],           n = "rlog") ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=names(DEGs_batch),                  MoreArgs = list(z=rld_df[,9:(8+nrow(samples))],                     n = "rlog") ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=paste0(names(DEGs_batch), "_mean"), MoreArgs = list(z=vsd_df_mean[,9:(9+length(conditions))],           n = "vst") ))
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=names(DEGs_batch),                  MoreArgs = list(z=vsd_df[,9:(8+nrow(samples))],                     n = "vst") ))
# shrunken log2FC (one value per condition already -- no separate "_mean" version to add).
# as.data.frame(): DEG_heatmap_ids's subset(z, ...) call assumes z is a
# data.frame (true of every other z argument here) -- on a bare matrix,
# subset() dispatches to subset.default (vector semantics), which would
# silently mangle row selection instead of properly picking rows.
invisible(mapply(FUN = DEG_heatmap_ids, x=DEGs_batch, y=names(DEGs_batch),                  MoreArgs = list(z=as.data.frame(log2FC_matrix),                     n = "log2FC") ))

# write DEG tables
DEG_write <- function(x, y) {
  write.table(subset(cts_summary_norm_mean, subset=row.names(cts_summary_norm_mean) %in% x), file=paste0(batch_dir, "batch_", y, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)
}
invisible(mapply(FUN = DEG_write, x=DEGs_batch, y=names(DEGs_batch)))
write.table(x=t(as.data.frame(lapply(DEGs_batch, FUN=length))), file=paste0(batch_dir, "DEGs_batch.tsv"), quote=F, sep="\t", row.names=T, col.names=F)

# ============================================================================
# One entry per feature (antisense run): sense and antisense calls of each
# comparison and of their union over comparisons, a feature called in both
# orientations counted once (one_entry_per_feature, antisense_calls.R).
# DEGs_both_orientations/<cond>_vs_<ref>_<class>.tsv, batch_<class>.tsv and
# DEGs_summary.tsv; Geneid carries _AS when the antisense orientation is kept.
# ============================================================================
if (strand == "_AS") {
  both_dir <- "DEGs_both_orientations/"
  dir.create(both_dir, showWarnings=FALSE)
  cond_means <- function(d) sapply(levels(d$condition), function(c) rowMeans(counts(d, normalized=TRUE)[, d$condition == c, drop=FALSE]))
  nm_AS <- cond_means(dds)
  rownames(nm_AS) <- paste0(rownames(nm_AS), "_AS")
  nm <- rbind(cond_means(dds_sense), nm_AS)
  res_both <- lapply(setNames(conditions, conditions), function(cond) {
    r_S <- results(dds_sense, contrast=c("condition", cond, args[2]))
    r_A <- all_res[[cond]]
    data.frame(Geneid=c(row.names(r_S), paste0(row.names(r_A), "_AS")), GeneId=c(row.names(r_S), row.names(r_A)),
               log2FoldChange=c(r_S$log2FoldChange, r_A$log2FoldChange), padj=c(r_S$padj, r_A$padj),
               called_up=c(is_call(r_S, SAF, "up"), all_kept[[cond]] & is_call(r_A, SAF, "up")),
               called_down=c(is_call(r_S, SAF, "down"), all_kept[[cond]] & is_call(r_A, SAF, "down")))
  })
  both_summary <- list()
  for (cl in c("TEGs", "TEs", "PCGs")) {
    ann_ids <- get(cl)$GeneId
    for (direction in c("up", "down")) {
      name <- paste0(direction, cl)
      calls <- lapply(res_both, function(r) r$Geneid[r[[paste0("called_", direction)]] & r$GeneId %in% ann_ids])
      changes <- lapply(setNames(conditions, conditions), orientation_change, norm_mean=nm, ref=args[2], direction=direction)
      for (cond in conditions) {
        e <- one_entry_per_feature(calls[[cond]], changes[[cond]])
        e <- merge(e, res_both[[cond]][, c("Geneid", "GeneId", "log2FoldChange", "padj")], by="Geneid", sort=FALSE)
        write.table(e, file=paste0(both_dir, cond, "_vs_", args[2], "_", name, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)
        both_summary[[name]][[cond]] <- nrow(e)
      }
      u <- union_one_entry_per_feature(calls, changes)
      u <- cbind(u, nm[u$Geneid, , drop=FALSE])
      write.table(u, file=paste0(both_dir, "batch_", name, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)
      both_summary[[name]][["batch"]] <- nrow(u)
    }
  }
  write.table(data.frame(DEGs=names(both_summary), do.call(rbind, lapply(both_summary, unlist)), check.names=FALSE),
              file=paste0(both_dir, "DEGs_summary.tsv"), quote=F, sep="\t", row.names=F, col.names=T)
}

unlink("Rplots.pdf") # delete that plot that is always created for some reason
