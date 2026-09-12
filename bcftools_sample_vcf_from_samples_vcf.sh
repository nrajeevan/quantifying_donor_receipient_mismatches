#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
known_sites=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf

#sample=ssssss1
sample=CFD-001_gDNA-Donor

gatk_dir=${work_dir}/results/gatk_new

# extract vcfs for a sample from vcfs of all samples.
EXTRACT_SAMPLE_SNPS_INDELS_FROM_ALL_SAMPLES=TRUE
if [ $EXTRACT_SAMPLE_SNPS_INDELS_FROM_ALL_SAMPLES ]; then
   all_samples_raw_snps=${gatk_dir}/genotypes_from_all_samples_GVCFs_AS_snps_recal.vcf.gz
   all_samples_raw_indels=${gatk_dir}/genotypes_from_all_samples_GVCFs_AS_indels_recal.vcf.gz
   raw_snps=${gatk_dir}/${sample}_raw_snps.vcf
   raw_indels=${gatk_dir}/${sample}_raw_indels.vcf

   bcftools view -s ${sample} ${all_samples_raw_snps} > ${raw_snps}
fi
