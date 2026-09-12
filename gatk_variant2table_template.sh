#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
known_sites=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf

sample=ssssss

results_dir=${work_dir}/results/gatk

if [ ! TRUE ];then
   data_dir=${work_dir}/results/bwa
   aligned_reads=${data_dir}/${sample}_sorted.bam
   sorted_dedup_reads=${results_dir}/${sample}_sorted_dedup_reads.bam
   recal_data_table=${results_dir}/${sample}_recal_data.table
   sorted_dedup_bqsr_reads=${results_dir}/${sample}_sorted_dedup_bqsr_reads.bam
   alignment_metrics=${results_dir}/${sample}_alignment_metrics.txt
   insert_size_metrics=${results_dir}/${sample}_insert_size_metrics.txt
   insert_size_histogram=${results_dir}/${sample}_insert_size_histogram.pdf
   raw_variants=${results_dir}/${sample}_raw_variants.vcf
   raw_snps=${results_dir}/${sample}_raw_snps.vcf
   raw_indels=${results_dir}/${sample}_raw_indels.vcf
fi

# table output
gatk_table_dir=${results_dir}/VCF_Tables

#raw_snps_table=${gatk_table_dir}/${sample}_raw_snps.tsv
#raw_indels_table=${gatk_table_dir}/${sample}_raw_indels.tsv



# Variant to Table method.
if [ TRUE ]; then

   if [ TRUE ];then
      declare -a vcf_files=(
         ${sample}_raw_snps.vcf ${sample}_raw_snps.tsv
         ${sample}_raw_indels.vcf ${sample}_raw_indels.tsv
         #${sample}_comm_snps.vcf ${sample}_comm_snps.tsv
         #${sample}_comm_indels.vcf ${sample}_comm_indels.tsv
      )

      for ((i=0;i< ${#vcf_files[@]} ;i+=2));
      do
         infile=${results_dir}/${vcf_files[i]}
         outfile=${results_dir}/VCF_Tables/${vcf_files[i+1]}
   
         gatk VariantsToTable \
            -R ${ref_fasta} \
            -V ${infile} \
            -F CHROM -F POS -F ID -F QUAL -F FILTER -F TYPE -F REF -F ALT -F AF -F AC -F DP -GF GT -GF DP -GF FT -GF GL -GF AD -GF AF -GF PL -GF GP -GF GQ -GF EC \
            --show-filtered \
            -O ${outfile}
      done
   fi

   if [ ! TRUE ];then
      declare -a vcf_files=(
         ${sample}-Donor_${sample}-Rec_snps_diff.vcf ${sample}_Donor-Rec_snps_diff.tsv
         ${sample}-Rec_${sample}-Donor_snps_diff.vcf ${sample}-Rec-Donor_snps_diff.tsv
         ${sample}-Donor_${sample}-Rec_indels_diff.vcf ${sample}-Donor-Rec_indels_diff.tsv
         ${sample}-Rec_${sample}-Donor_indels_diff.vcf ${sample}-Rec-Donor_indels_diff.tsv
      )
   fi

   if [ ! TRUE ];then
      declare -a vcf_files=(
         ${sample1}_snps_comm.vcf ${sample2}_snps_comm.tsv
         ${sample1}_indels_comm.vcf ${sample2}_indels_comm.tsv
      )
   fi

   if [ ! TRUE ];then
      for ((i=0;i< ${#vcf_files[@]} ;i+=2));
      do
         infile=${results_dir}/${vcf_files[i]}
         outfile=${results_dir}/VCF_Tables/${vcf_files[i+1]}
   
         gatk VariantsToTable \
            -R ${ref_fasta} \
            -V ${infile} \
            -F CHROM -F POS -F ID -F QUAL -F FILTER -F TYPE -F REF -F ALT -F AF -GF GT -GF DP -GF PL -GF GP -GF GQ \
            --show-filtered \
            -O ${outfile}
   
      done
   fi

#   gatk VariantsToTable \
#    -R ${ref_fasta} \
#    -V ${raw_snps} \
#    -F CHROM -F POS -F ID -F QUAL -F FILTER -F TYPE -F REF -F ALT -F AF -GF GT -GF DP -GF PL -GF GP -GF GQ \
#    --show-filtered \
#    -O ${raw_snps_table}
fi

