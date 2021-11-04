#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-00:10:00
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.out

# $1 is the folder with the raw counts (and a counts_summary file)
# $2 is the control condition prefix (eg Col0)
# $3 is an absolute path to the samples.txt table used for STARalign.sh script submission
# $4 is single-strand or paired-end
# $5 is quantseq or empty
# array size is the number of pairwise comparisons (number of conditions minus 1)
# NO HYPHEN, only special characters tolerated by R are . and _

# in the sample table, replicates must be indicated with R1, R2, R3...

sample=(`tail -n+2 $3 | cut -d "," -f2 | sort | uniq | grep -wv "$2" | sed -n ${SLURM_ARRAY_TASK_ID}p`)
echo -e "$(date) .. job ${SLURM_JOB_NAME} started on node $(/bin/hostname) with slurm ID: ${SLURM_JOB_ID} and array ${SLURM_ARRAY_TASK_ID}.\nfolder: ${1}\ncontrol: ${2}\nsample_table: ${3}\ntype: ${4} ${5}\nsample: ${sample}\n" 

ml build-env/f2020 r-bundle-bioconductor/3.9-foss-2019a-r-3.6.0

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    echo -e "\ndoing batch analysis... \n"
    Rscript DESeq2_from_featureCounts_batch_analysis.R $1 $2 $3 $4 sense $5 &
    Rscript DESeq2_from_featureCounts_batch_analysis.R $1 $2 $3 $4 antisense $5
    wait
fi

echo -e "\ndoing pairwise comparison... \n"

Rscript /groups/berger/user/pierre.bourguet/genomics/scripts/DESeq2/DESeq2_from_featureCounts_pairwise_comparison.R $1 $2 $sample $3 $4 $5

echo -e "\nplotting 100kb windows... \n"
module reset
ml build-env/2020 bedtools/2.27.1-foss-2018b r/4.0.2-foss-2018b
for i in ${1%/}/DEGs_${sample}_vs_${2}/${sample}_vs_${2}_up*mean.tsv ${1%/}/DEGs_${sample}_vs_${2}/${sample}_vs_${2}_down*mean.tsv ${1%/}/DEGs_${sample}_vs_${2}_AS/${sample}_vs_${2}_up*mean.tsv ${1%/}/DEGs_${sample}_vs_${2}_AS/${sample}_vs_${2}_down*mean.tsv ; do
	if [[ `wc -l $i | cut -f1 -d " "` -gt 1 ]]; then
        bedtools intersect -c -a /groups/berger/user/pierre.bourguet/genomics/Araport11/windows_chr1_5_100kb.tsv -b $i > ${i%mean.tsv}.tsv.100kbwindows_tmp # remove the "mean" in filename as it is useless
	    Rscript /groups/berger/user/pierre.bourguet/genomics/scripts/DESeq2/100kb_windows_plot.R ${i%mean.tsv}.tsv.100kbwindows_tmp
	    rm ${i%mean.tsv}.tsv.100kbwindows_tmp
    fi
done

echo "job finished at: $(date)"
