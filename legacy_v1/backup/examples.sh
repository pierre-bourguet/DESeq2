#!/usr/bin/env Rscript 

sbatch DESeq2_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/ Col0 ddm1_1st /groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_1st_7th_ribozero_STAR/samples.txt

sbatch DESeq2_from_kallisto_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/kallisto_output/kallistoData.Rdata WT cmt3 /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/ 3

sbatch DESeq2_from_kallisto_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/kallisto_output/kallistoData.Rdata WT cmt3 /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/ 3

sbatch DESeq2_from_kallisto_all_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/kallisto_output/kallistoData.Rdata WT /groups/berger/user/pierre.bourguet/genomics/scripts/bam_to_kallisto/info_h2aw2_ribovanish.tab /groups/berger/user/pierre.bourguet/genomics/RNAseq/h2aw-2_ribovanish/kallisto_SE/ 3
