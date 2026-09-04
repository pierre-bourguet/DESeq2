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
# $3 is an absolute path to the samples.txt table used for STAR submission

echo -e "job ${SLURM_JOB_NAME} started at: $(date) on node $(/bin/hostname) with slurm ID: ${SLURM_JOB_ID}." 

# you might need to change that depending on how replicates are named in the samples.txt table (rep1 or Rep1 or whatever)
samples=`cut -f1 $3 | sed 's/_rep.$//' | sort | uniq | grep -v "$2"`

for i in $samples; do
	sbatch /groups/berger/user/pierre.bourguet/genomics/scripts/DESeq2/DESeq2_from_featureCounts_pairwise_comparison.sbatch $1 $2 $i $3
done

echo "job finished at: $(date)"