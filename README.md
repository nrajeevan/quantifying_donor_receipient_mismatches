Files in this repository:

For any assistance in using any of these programs, please contact:<br>
	Nallakkandi Rajeevan, Ph.D. <br>
	Email: n.rajeevan@yale.edu

The following programs are written to ru in a linux cluster running a job scheduler.

First create a bash array of sample files as:
declare -a sample_array=(<br>
	Sample1_Donor_gdna<br>
	Sample1_Rec_gdna<br>
	Sample1_Rec_cfdna<br>
	Sample1_Rec_urdna<br>
	Sample2_Donor_gdna<br>
	Sample2_Rec_gdna<br>
	Sample2_Rec_cfdna<br>
	Sample2_Rec_urdna<br>
)

This array is used to create batch files and joblist for scheduling the jobs in the cluster.

1) FastQC quality checking

	FASTQC is aprogram to check the quality of the reads. This plots various quality metrics.
	
   a) fastqc_template.sh<br>
   b) make_fastqc_batch_files.sh<br>
   c) make_joblist_fastqc_from_template.sh<br>
   
2) BBDuk Adapter trim

	This program is used to adapter trim the raw R1 and R2 FASTQ files.
	
   a) bbduk_adapter_trim_template.sh<br>
   b) make_bbduk_batch_files.sh<br>
   c) make_joblist_bbduk_from_template.sh<br>

3) BWA alignment

	The BWA aligner is used to map paired-end reads to homo-sapien reference sequence.
	
   a) bwa_mem_align_template.sh<br>
   b) make_bwa_mem_align_batch_files.sh<br>
   c) make_joblist_bwa_mem_align.sh<br>

4) FastqScreen quality checking for contamination

   a) fastq_screen_template.sh<br>
   b) make_fastq_screen_batch_files.sh<br>
   c) make_joblist_fastq_screen_from_template.sh<br>

5) Variant call (diploid and quadraploid)

	We used GATK for variant call. Make sure the parameter -ploidy is set to 4 for tetraploid.
	
   a) gatk_variant_call_template.sh<br>
   b) make_gatk_variant_call_batch_files.sh<br>
   c) make_joblist_gatk_variant_call_from_template.sh<br>

6) Genotype GVCFs all samples

   a) gatk_combine_GVCFs_all_samples.sh<br>
   b) gatk_Genotype_GVCFs_all_samples.sh<br>
   c) make_joblist_gatk_Genotype_GVCFs_all_samples_from_template.sh<br>
   d) gatk_Filter_Variants_with_VQSR_all_samples.sh<br>
   e) bcftools_sample_vcf_from_samples_vcf.sh<br>
   
 7) Compute Donor/Recepient mismatch 
 
	These in-house developed programs are used to estimate various mismatch matrices.
	
   a) mismatch_in_cfDNA_or_urDNA.sh<br>
   b) mismatch_score_analysis_geneset.sh<br>
   c) mismatch_score_analysis_genome.sh<br>
   
