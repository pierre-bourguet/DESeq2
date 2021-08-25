#!/usr/bin/env bash

# this creates a summary file with number of DEGs for each genotype
# $1 is the path to the main output folder
cd ${1}
for i in `ls -d *DEGs`
    do
        sample=`echo $i | sed 's/_DEGs//'`
        if [ ! -f DEGs_summary.tsv ]; then
            echo -e "DEGs\t$sample" > DEGs_summary.tsv
            cat ${i}/${i}_summary.tsv >> DEGs_summary.tsv # create file if it does not exist yet
        else
            echo $sample > tmp
            cut -f2 ${i}/${i}_summary.tsv >> tmp
            paste DEGs_summary.tsv tmp > tmp2
            mv tmp2 DEGs_summary.tsv
            rm tmp
        fi
done