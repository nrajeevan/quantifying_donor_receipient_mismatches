#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
known_sites=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf

sample=ssssss

result_dir=${work_dir}/results/gatk
data_dir=${work_dir}/results/bwa
aligned_reads=${data_dir}/${sample}_sorted.bam
sorted_dedup_reads=${result_dir}/${sample}_sorted_dedup_reads.bam
recal_data_table=${result_dir}/${sample}_recal_data.table
sorted_dedup_bqsr_reads=${result_dir}/${sample}_sorted_dedup_bqsr_reads.bam
alignment_metrics=${result_dir}/${sample}_alignment_metrics.txt
insert_size_metrics=${result_dir}/${sample}_insert_size_metrics.txt
insert_size_histogram=${result_dir}/${sample}_insert_size_histogram.pdf
#dbsnp_file=${ref_dir}/GCF_000001405.40.gz
dbsnp_file=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf
raw_variants=${result_dir}/${sample}_raw_variants.vcf
raw_snps=${result_dir}/${sample}_raw_snps.vcf
raw_indels=${result_dir}/${sample}_raw_indels.vcf

if [ TRUE ]; then
   all_sample_snps_vcfs=${result_dir}/genotypes_from_all_samples_GVCFs_AS_snps_recal.vcf.gz
   all_sample_indels_vcfs=${result_dir}/genotypes_from_all_samples_GVCFs_AS_indels_recal.vcf.gz
   sample_snps_vcf=${result_dir}/${sample}_raw_snps.vcf
   sample_indels_vcf=${result_dir}/${sample}_raw_indels.vcf

   gatk SelectVariants -R ${ref_fasta} -V ${all_sample_snps_vcfs} --sample-name ${sample} -O ${sample_snps_vcf}
   gatk SelectVariants -R ${ref_fasta} -V ${all_sample_indels_vcfs} --sample-name ${sample} -O ${sample_indels_vcf}
fi
