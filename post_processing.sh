#!/usr/bin/env bash
# Manifest-driven replacement for the old post_processing.sh.
#
# Builds DEGs_summary.tsv (sense) and DEGs_AS_summary.tsv (antisense) from
# the DEGs_manifest[.tsv|_AS.tsv] files written by DESeq2_pipeline.R, instead
# of globbing "DEGs_*" in the output directory. This means a stale DEGs_*
# folder left over from an earlier run (different thresholds, a removed
# condition, a leftover test...) can never silently corrupt the summary: only
# folders actually listed in the current manifest are read, and the summary
# file is rebuilt from scratch every time (safe to rerun anytime).
#
# usage: post_processing.sh <output_dir>
# <output_dir> is the DESeq2 output directory (args[7] passed to
# DESeq2_pipeline.R, or args[1] if that was omitted).

set -euo pipefail

output_dir=$1
cd "$output_dir"

build_summary () { # $1 = manifest file, $2 = summary file to (re)write, $3 = strand suffix ("" or "_AS")
  local manifest=$1 summary=$2 suffix=$3
  if [[ ! -s "$manifest" ]]; then
    echo "post_processing.sh: no ${manifest} (or empty) in ${output_dir}, skipping ${summary}"
    return 0
  fi

  # header row: "DEGs" followed by one column per comparison in the manifest
  # (dir names carry the strand suffix, eg "DEGs_cond_vs_ref_AS"; the files
  # inside a comparison's folder never carry it -- only the folder does)
  { printf "DEGs"; while IFS= read -r dir; do stem="${dir#DEGs_}"; printf "\t%s" "${stem%$suffix}"; done < "$manifest"; printf "\n"; } > "$summary"

  # one row per DEG category, one value per comparison
  for category in upTEGs upTEs upPCGs downTEGs downTEs downPCGs; do
    row="$category"
    while IFS= read -r dir; do
      stem="${dir#DEGs_}" ; stem="${stem%$suffix}"
      value=$(awk -F'\t' -v c="$category" '$1==c{print $2}' "${dir}/${stem}_DEGs_summary.tsv")
      row="${row}"$'\t'"${value:-NA}"
    done < "$manifest"
    printf "%s\n" "$row" >> "$summary"
  done
}

build_summary "DEGs_manifest.tsv"    "DEGs_summary.tsv"    ""
build_summary "DEGs_manifest_AS.tsv" "DEGs_AS_summary.tsv" "_AS"

# cross-comparison DEG intersections (sense and antisense, when their manifest exists)
"$(dirname "$0")/intersect_DEG_batch.sh" "$output_dir" ""
"$(dirname "$0")/intersect_DEG_batch.sh" "$output_dir" "_AS"
