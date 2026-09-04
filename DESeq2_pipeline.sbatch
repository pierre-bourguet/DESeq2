#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-03:00:00
#SBATCH --output=/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/logs/%A_%x.out
#SBATCH --error=/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2/logs/%A_%x.out

# $1 is the mapping-pipeline folder (read-only input: counts_summary.tsv / counts_AS_summary.tsv)
# $2 is the control/reference condition prefix (eg Col0, WT)
# $3 is an absolute path to the samples.txt table used for STARalign.sbatch submission
# $4 is single-strand or paired-end
# $5 is quantseq or empty
# $6 is optional: a dedicated output directory for all DESeq2 outputs.
#    Leave empty (or omit) to keep the old behavior of writing into $1.
#    Recommended for new runs: point $1 at the untouched mapping folder and
#    $6 at a fresh results folder -- no more `cp`-ing mapping outputs around
#    before (re-)running DESeq2, since $1 is never written to.
#
# This is no longer a --array job: DESeq2_pipeline.R fits ONE Master model
# per strand and internally loops over every non-reference condition to
# extract all pairwise comparisons, so a single sense run + a single
# antisense run cover the whole experiment. Post-processing (DEG summary +
# cross-condition intersects) now runs automatically at the end, reading the
# manifest files DESeq2_pipeline.R writes -- no separate manual step needed.
#
# NO HYPHEN in condition/sample names, only "." and "_" are tolerated by R.
# In the sample table, replicates must be indicated with R1, R2, R3...

# Hardcoded, not derived from $0: sbatch copies the submitted script into a
# spool dir (/var/spool/slurm/job<id>/slurm_script) and runs it from there,
# so "$(dirname "$0")" resolves to the spool dir at runtime, not to where
# this script actually lives -- DESeq2_pipeline.R / post_processing.sh would
# silently fail to be found otherwise.
pipeline_dir="/groups/berger/user/pierre.bourguet/shared/pipelines/DESeq2"
inputdir=${1%/}/ # deal with cases where the user did not include a trailing slash in $1
outputdir=${6:-$1} ; outputdir=${outputdir%/}/
echo -e "$(date) .. job ${SLURM_JOB_NAME} started on node $(/bin/hostname) with slurm ID: ${SLURM_JOB_ID}.\ninput (mapping) folder: ${1}\noutput (DESeq2) folder: ${outputdir}\ncontrol: ${2}\nsample_table: ${3}\ntype: ${4} ${5}\n"

module load build-env/f2022 r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1

# Without this, the R module's linked BLAS/OpenMP libraries spawn as many
# threads per process as the node has cores, ignoring --cpus-per-task=1 --
# with 2 concurrent R processes (sense + antisense) that means severe core
# oversubscription/contention rather than proportional speedup. Pin every
# threading knob to 1 so each process actually uses only its allocated CPU.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

echo -e "\nrunning DESeq2 Master model + pairwise extraction (sense & antisense in parallel)...\n"
Rscript "${pipeline_dir}/DESeq2_pipeline.R" "$inputdir" "$2" "$3" "$4" sense "$5" "$6" &
Rscript "${pipeline_dir}/DESeq2_pipeline.R" "$inputdir" "$2" "$3" "$4" antisense "$5" "$6"
wait

echo -e "\nrunning post-processing (DEGs_summary + cross-condition intersects, sense & antisense)...\n"
module reset
ml build-env/2020 r/4.0.2-foss-2018b
"${pipeline_dir}/post_processing.sh" "$outputdir"

echo -e "\nplotting 100kb windows for every comparison (sense & antisense)...\n"
ml bedtools/2.27.1-foss-2018b

# manifests list exactly the DEGs_<cond>_vs_<ref>[_AS] folders produced this
# run -- read from there instead of re-deriving them from the sample table,
# so this loop can never drift from what DESeq2_pipeline.R actually created.
for manifest in "${outputdir}DEGs_manifest.tsv" "${outputdir}DEGs_manifest_AS.tsv"; do
  [[ -s "$manifest" ]] || continue
  while IFS= read -r comparison; do
    dir="${outputdir}${comparison}"
    stem="${comparison#DEGs_}" ; [[ "$manifest" == *_AS.tsv ]] && stem="${stem%_AS}"
    for i in "$dir"/${stem}_up*mean.tsv "$dir"/${stem}_down*mean.tsv; do
      [[ -e "$i" ]] || continue
      if [[ $(wc -l < "$i") -gt 1 ]]; then
        bedtools intersect -c -a /groups/berger/user/pierre.bourguet/genomics/Araport11/windows_chr1_5_100kb.tsv -b "$i" > "${i%mean.tsv}.tsv.100kbwindows_tmp" # drop "mean" from filename, it is no longer useful
        Rscript "${pipeline_dir}/100kb_windows_plot.R" "${i%mean.tsv}.tsv.100kbwindows_tmp"
        rm "${i%mean.tsv}.tsv.100kbwindows_tmp"
      fi
    done
  done < "$manifest"
done

echo "job finished at: $(date)"
