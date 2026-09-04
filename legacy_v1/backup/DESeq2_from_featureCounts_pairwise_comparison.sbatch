#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=0-02:00:00
#SBATCH --output=logs/%A_%x.out
#SBATCH --error=logs/%A_%x_error.out

# $1 is the folder with the raw counts (and a counts_summary file)
# $2 is the control condition prefix (eg Col0)
# $3 is the treatment condition prefix
# $4 is an absolute path to the samples.txt table used for STAR submission

ml r/4.0.2-foss-2018b build-env/f2020 r-bundle-bioconductor/3.9-foss-2019a-r-3.6.0

Rscript /groups/berger/user/pierre.bourguet/genomics/scripts/DESeq2/DESeq2_from_featureCounts_pairwise_comparison.R $1 $2 $3 $4