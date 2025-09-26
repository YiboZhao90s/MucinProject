# Use vcftools to calculate MAFs and generate genotype matrix
module load vcftools/0.1.14-uowdqmf  
vcftools --gzvcf merged_MUC_sorted_dedup_MUC5AC.vcf.gz --freq --out allele_freq_MUC5AC
vcftools --gzvcf merged_MUC_sorted_dedup_MUC5AC.vcf.gz --012 --out genotype_matrix_MUC5AC
# Note: for genotype matrix, rows are samples and variant are columns
