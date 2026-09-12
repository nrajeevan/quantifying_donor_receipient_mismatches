Files in this repository:


1. BBDuk Adapter trim

   a) bbduk_adapter_trim_template.sh<br>
   b) make_bbduk_batch_files.sh
   c) make_joblist_bbduk_from_template.sh

2) BWA alignment

   a) bwa_mem_align_template.sh
   b) make_bwa_mem_align_batch_files.sh
   c) make_joblist_bwa_mem_align.sh

3) FastQC quality checking

   a) fastqc_template.sh
   b) make_fastqc_batch_files.sh
   c) make_joblist_fastqc_from_template.sh

4) FastqScreen quality checking for contamination

   a) fastq_screen_template.sh
   b) make_fastq_screen_batch_files.sh
   c) make_joblist_fastq_screen_from_template.sh

5) Variant call (diploid and quadraploid)

   a) gatk_variant_call_template.sh
   b) make_gatk_variant_call_batch_files.sh
   c) make_joblist_gatk_variant_call_from_template.sh

6) Genotype GVCFs all samples

   a) gatk_combine_GVCFs_all_samples.sh
   b) gatk_Genotype_GVCFs_all_samples.sh
   c) make_joblist_gatk_Genotype_GVCFs_all_samples_from_template.sh
   d) gatk_Filter_Variants_with_VQSR_all_samples.sh
   e) bcftools_sample_vcf_from_samples_vcf.sh
   
 7) Compute Donor/Recepient mismatch 
 
   a) mismatch_in_cfDNA_or_urDNA.sh
   b) mismatch_score_analysis_geneset.sh
   c) mismatch_score_analysis_genome.sh
   
