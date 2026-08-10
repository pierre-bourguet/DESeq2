sbatch --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end

sbatch -p m --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end

sbatch -p m --array=1-8 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/rerun/w2_h1_suvh456_cmt3_ribovanish/ WT /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w2_h1_suvh456_cmt3_ribovanish.samples paired-end

sbatch -p m --array=1-3 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/2020_Rougee_ddm1_clf/ Col /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/2020_Rougee_ddm1_clf.samples single-strand quantseq

sbatch -p m --array=1-3 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/w_heat_stress_qseq/ WT_23 /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/w_heat_stress_qseq.samples single-strand quantseq

for i in 1M 2M 4M 500K ; do sbatch -p g --array=1-5 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/ddm1_hetero_kanno_tagseq_in_vitro_whole_seedlings_10_14d_coverage_testing/500K_reads Col_10 /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/ddm1_kanno_500K_sample.samples single-strand quantseq ; done

sbatch --array=1-18 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/tagseq_03_cdca7 WT /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/sample_lists/DESeq2_tagseq_03_cdca7.tsv single-strand quantseq

sbatch --array=1-3 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/genomics/RNAseq/MM2D_aphidicolin_x5_PE75_no_ctrl Aph_0h /groups/berger/user/pierre.bourguet/genomics/RNAseq/sample_tables/MM2D_aphidicolin_Anna_no_ctrl.samples paired-end

sbatch -p m --array=1-5 DESeq2_from_featureCounts_pairwise_comparison.sbatch /groups/berger/user/pierre.bourguet/projects/03_regulatory_TEs/results/syringolinA_Bonnet_2023/DESeq2/ Col_ctrl /groups/berger/user/pierre.bourguet/projects/03_regulatory_TEs/data/syringolinA_Bonnet_2023/DESeq2_formatting/GSE174350_RNAseq_counts.samples paired-end
