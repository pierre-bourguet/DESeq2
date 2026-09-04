#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-00:30:00
#SBATCH --output=/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/logs/%A_%x.out
#SBATCH --error=/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/logs/%A_%x.out

# Wrapper for DESeq2_pairwise_contrasts.R -- extracts EXTRA pairwise contrasts
# (e.g. combined double-mutant vs single-mutant) from an ALREADY-FITTED Master
# dds produced by a prior DESeq2_pipeline.sbatch run, without re-estimating
# dispersion (see DESeq2_pairwise_contrasts.R header for the rationale). Much
# lighter than DESeq2_pipeline.sbatch -- no dispersion fit, no PCA/rlog/vst/MoR
# batch exports (those already exist from the original run) -- hence the
# smaller mem/time request.
#
# $1 = mapping-pipeline folder (read-only input: counts_summary.tsv / counts_AS_summary.tsv)
# $2 = absolute path to the sample table used for the ORIGINAL master run
# $3 = quantseq or empty (must match the normalization of the original run)
# $4 = dds_dir: READ-ONLY output directory of the completed DESeq2_pipeline.R
#      run for this dataset (holds dds[_AS].rds). Never written to.
# $5 = comma-separated combined/double-mutant conditions, no spaces
# $6 = comma-separated single-mutant baselines, same length/order as $5
# $7 = optional: dedicated output directory for the new
#      DEGs_<combined>_vs_<single>[_AS]/ folders + DEGs_manifest_custom[_AS].tsv.
#      Leave empty (or omit) to write alongside the original run in $4.
#
# NO HYPHEN in condition/sample names, only "." and "_" are tolerated by R.

pipeline_dir="/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2"
inputdir=${1%/}/
ddsdir=${4%/}/
outputdir=${7:-$4} ; outputdir=${outputdir%/}/
echo -e "$(date) .. job ${SLURM_JOB_NAME} started on node $(/bin/hostname) with slurm ID: ${SLURM_JOB_ID}.\ninput (mapping) folder: ${1}\ndds (read-only) folder: ${ddsdir}\noutput folder: ${outputdir}\nsample_table: ${2}\ncombined: ${5}\nsingle (baseline): ${6}\n"

module load build-env/f2022 r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

echo -e "\nrunning DESeq2 custom pairwise contrasts (sense & antisense in parallel)...\n"
Rscript "${pipeline_dir}/DESeq2_pairwise_contrasts.R" "$inputdir" "$2" sense "$3" "$ddsdir" "$5" "$6" "$outputdir" &
Rscript "${pipeline_dir}/DESeq2_pairwise_contrasts.R" "$inputdir" "$2" antisense "$3" "$ddsdir" "$5" "$6" "$outputdir"
wait

echo -e "\nbuilding DEGs_custom_summary.tsv from the manifest(s)...\n"
cd "$outputdir"
build_summary () { # $1 = manifest file, $2 = summary file to (re)write, $3 = strand suffix ("" or "_AS")
  local manifest=$1 summary=$2 suffix=$3
  if [[ ! -s "$manifest" ]]; then
    echo "no ${manifest} (or empty) in ${outputdir}, skipping ${summary}"
    return 0
  fi
  { printf "DEGs"; while IFS= read -r dir; do stem="${dir#DEGs_}"; printf "\t%s" "${stem%$suffix}"; done < "$manifest"; printf "\n"; } > "$summary"
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
build_summary "DEGs_manifest_custom.tsv"    "DEGs_custom_summary.tsv"    ""
build_summary "DEGs_manifest_custom_AS.tsv" "DEGs_custom_AS_summary.tsv" "_AS"

echo "job finished at: $(date)"
