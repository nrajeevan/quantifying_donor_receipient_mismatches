#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
known_sites=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf
dbsnp_file=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf
#dbsnp_file=${ref_dir}/GCF_000001405.40.gz

sample=ssssss
chrom=cccccc

# original sample R1 R2 files
#data_dir=${home_dir}/Data/Sample_${sample}

# bbduk trimmed sample R1 R2 files
# data_dir=${work_dir}/results/bbduk

results_dir=${work_dir}/results/gatk
data_dir=${work_dir}/results/bwa
aligned_reads=${data_dir}/${sample}_aligned_sorted_RG.bam
sorted_dedup_reads=${results_dir}/${sample}_sorted_dedup_reads.bam
sorted_dedup_bqsr_reads=${results_dir}/${sample}_sorted_dedup_bqsr_reads.bam
recal_data_table=${results_dir}/${sample}_recal_data.table
alignment_metrics=${results_dir}/${sample}_alignment_metrics.txt
insert_size_metrics=${results_dir}/${sample}_insert_size_metrics.txt
insert_size_histogram=${results_dir}/${sample}_insert_size_histogram.pdf

# without allele specific annotation
#raw_gvcf_variants=${results_dir}/${sample}_raw_variants.g.vcf.gz

# with allele specific annotations
raw_gvcf_variants=${results_dir}/${sample}_raw_variants_AS.g.vcf.gz

raw_variants=${results_dir}/${sample}_raw_variants.vcf
raw_snps=${results_dir}/${sample}_raw_snps.vcf
raw_indels=${results_dir}/${sample}_raw_indels.vcf

if [ TRUE ]; then
   
   # temporarily blocking on 6/9/2025
   if [ ! TRUE ]; then
      # Mark duplicates
      # temporary unmask this 
      gatk MarkDuplicatesSpark -I ${aligned_reads} -O ${sorted_dedup_reads}
   
      # Mark Duplicates and Sort
      gatk BaseRecalibrator -I ${sorted_dedup_reads} -R ${ref_fasta} --known-sites ${known_sites} -O ${recal_data_table}
      gatk ApplyBQSR -I ${sorted_dedup_reads} -R ${ref_fasta} --bqsr-recal-file ${recal_data_table} -O ${sorted_dedup_bqsr_reads}
   
      echo "STEP 5: Collect Alignment & Insert Size Metrics"
   
      gatk CollectAlignmentSummaryMetrics R=${ref_fasta} I=${sorted_dedup_bqsr_reads} O=${alignment_metrics}
      gatk CollectInsertSizeMetrics INPUT=${sorted_dedup_bqsr_reads} OUTPUT=${insert_size_metrics} HISTOGRAM_FILE=${insert_size_histogram}
   fi

   # for single sample - don't use this, instead, use the following to generate g.vcf output.
   if [ ! TRUE ]; then
      echo "STEP 6: Call Variants - gatk haplotype caller"

      gatk HaplotypeCaller -R ${ref_fasta} --dbsnp ${dbsnp_file} -I ${sorted_dedup_bqsr_reads} -O ${raw_variants}

      # extract SNPs & INDELS

      gatk SelectVariants -R ${ref_fasta} -V ${raw_variants} --select-type SNP -O ${raw_snps}
      gatk SelectVariants -R ${ref_fasta} -V ${raw_variants} --select-type INDEL -O ${raw_indels}
   fi

   # for multi-sample - generates g.vcf output for each sample and then use GenotypeGVCFs.
   if [ TRUE ]; then
      echo "STEP 6: Call Variants - gatk haplotype caller"

      # without allele speific annaotation
      if [ ! TRUE ]; then
         gatk HaplotypeCaller \
            -R ${ref_fasta} \
	    --dbsnp ${dbsnp_file} \
	    -I ${sorted_dedup_bqsr_reads} \
	    -O ${raw_gvcf_variants} \
	    -ERC GVCF
      fi

      # with allele speific annaotation
      # Reanalyzing with Ploidy = 4.
      if [ TRUE ]; then
         gatk HaplotypeCaller \
            -R ${ref_fasta} \
	    --dbsnp ${dbsnp_file} \
	    -I ${sorted_dedup_bqsr_reads} \
	    -O ${raw_gvcf_variants} \
	    -ERC GVCF \
            --native-pair-hmm-threads 16 \
	    -G StandardAnnotation \
	    -G AS_StandardAnnotation
      fi
   fi
fi

