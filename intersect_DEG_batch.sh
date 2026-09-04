#!/usr/bin/env bash
# Manifest-driven replacement for the old intersect_DEG_batch.sh.
#
# For every pair of comparisons listed in the manifest, counts the overlap
# (by Geneid, column 7) between their DEG lists, for each of the 6 DEG
# categories, and renders a heatmap per comparison + one heatmap per DEG
# category summarizing all comparisons pairwise. Reads DEGs_manifest[.tsv|
# _AS.tsv] instead of globbing "DEGs_*", and uses a private mktemp -d
# workspace instead of bare tmp/header/tmp_j/tmp_k files in the working
# directory -- the old script's plain filenames were a real collision risk
# now that sense and antisense DESeq2_pipeline.R runs (and therefore two
# post_processing.sh invocations) can happen back-to-back in the same folder.
#
# usage: intersect_DEG_batch.sh <output_dir> [strand_suffix]
#   strand_suffix: "" for sense (default) or "_AS" for antisense

set -euo pipefail

output_dir=$1
suffix=${2:-}
cd "$output_dir"

manifest="DEGs_manifest${suffix}.tsv"
if [[ ! -s "$manifest" ]]; then
  echo "intersect_DEG_batch.sh: no ${manifest} (or empty) in ${output_dir}, skipping intersect${suffix}"
  exit 0
fi

ml build-env/2020 r/4.0.2-foss-2018b

outdir="intersect_DEG${suffix}"
mkdir -p "${outdir}/tables"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

mapfile -t dirs < "$manifest"

intersect () { # $1, $2: two DEG tsv files; counts shared Geneids (column 7)
  comm -12 <(tail -n+2 "$1" | cut -f7 | sort) <(tail -n+2 "$2" | cut -f7 | sort) | wc -l
}

# short label per comparison (mutant name only, "_vs_ref[_AS]" stripped) for headers
touch "$work/header"
for dir in "${dirs[@]}"; do
  stem="${dir#DEGs_}" ; stem="${stem%$suffix}"
  label="${stem%_vs_*}"
  echo "$label" | paste "$work/header" - > "$work/tmp" && mv "$work/tmp" "$work/header"
done

# one heatmap table per comparison: how much does it overlap every other comparison, per DEG category
for dir_i in "${dirs[@]}"; do
  stem_i="${dir_i#DEGs_}" ; stem_i="${stem_i%$suffix}"
  : > "$work/tmp_j" # truncate/create -- NOT touch: tmp_j must start empty each outer iteration, and touch does not clear an existing file
  for dir_j in "${dirs[@]}"; do
    stem_j="${dir_j#DEGs_}" ; stem_j="${stem_j%$suffix}"
    : > "$work/tmp_k" # truncate/create -- same reason as tmp_j above
    for k in _upTEGs.tsv _upTEs.tsv _upPCGs.tsv _downTEGs.tsv _downTEs.tsv _downPCGs.tsv; do
      ref_file="${dir_i}/${stem_i}${k}"
      file="${dir_j}/${stem_j}${k}"
      intersect "$ref_file" "$file" | cat "$work/tmp_k" - > "$work/tmp" && mv "$work/tmp" "$work/tmp_k"
    done
    paste "$work/tmp_j" "$work/tmp_k" > "$work/tmp" && mv "$work/tmp" "$work/tmp_j"
  done
  cat "$work/header" "$work/tmp_j" | sed 's/\t//' > "$work/tmp_i"
  echo -e "DEG\nupTEG\nupTE\nupPCG\ndownTEG\ndownTE\ndownPCG" | paste - "$work/tmp_i" > "${outdir}/${stem_i}_intersect_with_others.tsv"
  Rscript /groups/berger/user/pierre.bourguet/genomics/scripts/intersect_DEG/table_heatmap.R "${outdir}/${stem_i}_intersect_with_others.tsv"
  mv "${outdir}/${stem_i}_intersect_with_others.tsv" "${outdir}/tables/"
  [[ -f "${outdir}/${stem_i}_intersect_with_others.pdf" ]] && mv "${outdir}/${stem_i}_intersect_with_others.pdf" "${outdir}/tables/"
done

# one heatmap table per DEG category: all comparisons vs all comparisons
for k in upTEGs upTEs upPCGs downTEGs downTEs downPCGs; do
  : > "$work/tmp_j" # truncate/create -- NOT touch: tmp_j must start empty each outer iteration, and touch does not clear an existing file
  for dir_i in "${dirs[@]}"; do
    stem_i="${dir_i#DEGs_}" ; stem_i="${stem_i%$suffix}"
    : > "$work/tmp_k" # truncate/create -- same reason as tmp_j above
    for dir_j in "${dirs[@]}"; do
      stem_j="${dir_j#DEGs_}" ; stem_j="${stem_j%$suffix}"
      ref_file="${dir_i}/${stem_i}_${k}.tsv"
      file="${dir_j}/${stem_j}_${k}.tsv"
      intersect "$ref_file" "$file" | paste "$work/tmp_k" - > "$work/tmp" && mv "$work/tmp" "$work/tmp_k"
    done
    cat "$work/tmp_j" "$work/tmp_k" > "$work/tmp" && mv "$work/tmp" "$work/tmp_j"
  done
  cat "$work/header" "$work/tmp_j" | sed 's/\t//' > "$work/tmp_i"
  tr "\t" "\n" < "$work/header" | paste - "$work/tmp_i" > "${outdir}/${k}_all_intersect.tsv"
  Rscript /groups/berger/user/pierre.bourguet/genomics/scripts/intersect_DEG/table_heatmap.R "${outdir}/${k}_all_intersect.tsv"
  mv "${outdir}/${k}_all_intersect.tsv" "${outdir}/tables/"
  [[ -f "${outdir}/${k}_all_intersect.pdf" ]] && mv "${outdir}/${k}_all_intersect.pdf" "${outdir}/tables/"
done
