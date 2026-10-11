# Antisense DEG call rules, shared by DESeq2_pipeline.R, DESeq2_pairwise_contrasts.R and the 3-tag-rna-seq pipeline (deseq2/DESeq2_tagseq.R).
#
# An antisense feature counts the reads on the strand opposite to its annotation. It is fitted apart from the sense features (own dispersion trend,
# own BH adjustment) with the size factors of the sense fit, on raw counts. An antisense up or down call is kept only if it passes both rules:
#   leakage                  in the libraries where antisense is higher (the compared condition for an up call, the reference for a down call),
#                            the normalized antisense count exceeds leak_fold (10) times the strand leakage expected from the sense count of the
#                            same feature (library leak rate x sense count). Leak rate of a library: median antisense/sense count ratio at the
#                            2,000 protein-coding genes with the highest mean normalized sense count.
#   opposite-strand partner  no feature overlapping it (>= 1 bp) on the opposite strand is called in sense, in the same direction, in the same
#                            contrast (the antisense reads would be that feature's transcript).
# Full results tables keep every feature with the two values the rules use: antisense_over_leakage and opposite_strand_sense_DEG.
# Lists built from the calls count a feature called in both orientations once ("sense wins", per condition or across conditions):
# pipelines/3-tag-rna-seq/README.md, "Sense and antisense features".
# Benchmark and simulation behind the rules: shared/tools/deseq2_antisense_benchmark/ (results in the projects listed in tools_inventory.tsv).

leak_fold <- 10

# rows `ids` of matrix M, zero-filled where M has no such row
rows_or_zero <- function(M, ids) {
  out <- matrix(0, length(ids), ncol(M), dimnames = list(ids, colnames(M)))
  k <- ids %in% rownames(M)
  out[k, ] <- M[ids[k], , drop = FALSE]
  out
}

# S, A: raw sense and antisense counts (features x libraries, antisense rows named by the sense ID); sf: sense size factors; pcg: protein-coding IDs
leak_rates <- function(S, A, sf, pcg, n_top = 2000) {
  pcg <- intersect(pcg, rownames(S))
  top <- head(pcg[order(-rowMeans(sweep(S[pcg, , drop = FALSE], 2, sf[colnames(S)], "/")))], n_top)
  apply(rows_or_zero(A, top)[, colnames(S), drop = FALSE] / pmax(S[top, , drop = FALSE], 1), 2, median)
}

# Normalized antisense count over expected leakage, summed over the libraries of the condition (log2FC >= 0) or of the reference (log2FC < 0).
# Inf when no leakage is expected, NA when the feature has no antisense count in those libraries.
antisense_over_leakage <- function(S, A, sf, leak, ids, lfc, cond_libs, ref_libs) {
  ratio <- function(libs) {
    num <- rowSums(sweep(rows_or_zero(A, ids)[, libs, drop = FALSE], 2, sf[libs], "/"))
    den <- rowSums(sweep(rows_or_zero(S, ids)[, libs, drop = FALSE], 2, leak[libs] / sf[libs], "*"))
    ifelse(num == 0, NA_real_, num / den)
  }
  unname(ifelse(is.na(lfc) | lfc >= 0, ratio(cond_libs), ratio(ref_libs)))
}

# Pairs of features overlapping (>= 1 bp) on opposite strands. ann: data.frame with id, chr, start, end, strand
opposite_strand_pairs <- function(ann) {
  gr <- GenomicRanges::GRanges(ann$chr, IRanges::IRanges(ann$start, ann$end), strand = ann$strand)
  hits <- GenomicRanges::findOverlaps(gr, gr, ignore.strand = TRUE)
  q <- S4Vectors::queryHits(hits); s <- S4Vectors::subjectHits(hits)
  opposite <- as.character(GenomicRanges::strand(gr))[q] != as.character(GenomicRanges::strand(gr))[s]
  data.frame(id = ann$id[q[opposite]], partner = ann$id[s[opposite]])
}

# Opposite-strand partners called in sense in the direction of each antisense log2FC: comma-separated IDs, NA if none
opposite_strand_sense_calls <- function(ids, lfc, pairs, sense_up, sense_down) {
  partners <- function(called) {
    p <- pairs[pairs$partner %in% called, ]
    vapply(split(p$partner, p$id), paste, character(1), collapse = ",")
  }
  up <- partners(sense_up); down <- partners(sense_down)
  unname(ifelse(is.na(lfc) | lfc >= 0, up[ids], down[ids]))
}

antisense_call_kept <- function(over_leakage, partners) !is.na(over_leakage) & over_leakage > leak_fold & is.na(partners)
