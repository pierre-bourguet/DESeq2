#!/usr/bin/env bash

# this creates a summary file with number of DEGs for each genotype
# $1 is the path to the main output folder
current_dir=$PWD
cd ${1}
for i in `ls -d DEGs_* | grep -v "_AS" | grep -v "_batch"`
    do
        sample=`echo $i | sed 's/DEGs_//'`
        if [ ! -f DEGs_summary.tsv ]; then
            echo -e "DEGs\t$sample" > DEGs_summary.tsv
            cat ${i}/${sample}_DEGs_summary.tsv >> DEGs_summary.tsv # create file if it does not exist yet
        else
            echo $sample > tmp
            cut -f2 ${i}/${sample}_DEGs_summary.tsv >> tmp
            paste DEGs_summary.tsv tmp > tmp2
            mv tmp2 DEGs_summary.tsv
            rm tmp
        fi
done

for i in `ls -d DEGs_*_AS | grep -v "_batch"`
    do
        sample=`echo $i | sed 's/DEGs_//'`
        if [ ! -f DEGs_AS_summary.tsv ]; then
            echo -e "DEGs\t$sample" > DEGs_AS_summary.tsv
            cat ${i}/${sample}_DEGs_summary.tsv >> DEGs_AS_summary.tsv # create file if it does not exist yet
        else
            echo $sample > tmp
            cut -f2 ${i}/${sample}_DEGs_summary.tsv >> tmp
            paste DEGs_AS_summary.tsv tmp > tmp2
            mv tmp2 DEGs_AS_summary.tsv
            rm tmp
        fi
done

# generates different files to intersect DEGs dectected in samples
bash /groups/berger/user/pierre.bourguet/genomics/scripts/intersect_DEG/intersect_DEG_batch.sh

cd $current_dir