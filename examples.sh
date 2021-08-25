#!/usr/bin/env Rscript 

sbatch --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end

sbatch -p m --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end

sbatch -p m --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end
