#!/usr/bin/env Rscript
# ============================================================================
# install_apeglm.R
#
# One-off, MANUAL install of the apeglm Bioconductor package into a
# project-local library, since apeglm is not bundled in any
# r-bundle-bioconductor module on this cluster (checked 3.8 through 3.19,
# the newest available -- see the comment near the top of
# DESeq2_pipeline.R).
#
# Run this ONCE, manually, from RStudio (or an interactive R session on a
# compute node -- never from a login-node shell, never via Rscript from
# outside RStudio). Load the same module the pipeline itself uses first, so
# the R version and compiler toolchain match:
#
#   module load build-env/f2022 r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1
#
# apeglm ships C++ code (Rcpp/RcppEigen) that needs to compile at install
# time; the gcc toolchain that comes with the module above should be enough.
#
# Installs into a project-local library folder (Rlib/, next to this script)
# instead of your personal ~/R library, so it's shared with anyone using
# this pipeline and never touches the central module.
# ============================================================================

lib_dir <- "/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/Rlib"
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

# BiocManager itself ships inside r-bundle-bioconductor already, so this
# should just find it -- the install.packages() fallback only fires if it's
# somehow missing.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib_dir, repos = "https://cloud.r-project.org")
}

BiocManager::install(
  "apeglm",
  lib    = lib_dir,
  update = FALSE,
  ask    = FALSE
)

cat(
  "\nDone. To use apeglm from any script afterwards, add this line BEFORE\n",
  "library(apeglm) (or before lfcShrink(..., type=\"apeglm\")):\n\n",
  "  .libPaths(c(\"", lib_dir, "\", .libPaths()))\n\n",
  sep = ""
)
