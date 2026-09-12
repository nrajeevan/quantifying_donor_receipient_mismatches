#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

. ${work_dir}/sample_arrays.sh


# CombineGVCFs where each sample GVCFs were generated WITHOUT allele speific annaotation
if [ ! TRUE ]; then
   samples=(${samples_orig[@]})

   combined_variants=""
   for sample in "${samples[@]}"
   do
   #	echo "--variant ${result_dir}/${sample}_raw_variants.g.vcf.gz \ "
   	combined_variants+="--variant ${result_dir}/${sample}_raw_variants.g.vcf.gz "
   done
   
   #echo ${!combined_variants}
   
   gatk_command="gatk CombineGVCFs \
      -R ${ref_fasta} \
      ${combined_variants} \
      -O ${result_dir}/combined_samples_raw_variants.g.vcf.gz"
   
   echo "gatk_command: ${gatk_command}"
   
   eval "${gatk_command}"
fi   

# CombineGVCFs where each sample GVCFs were generated WITH allele speific annaotation
# This section for control samples: samples_C1_C2 in sample_arrays.sh
# declare -a samples_orig=(
# )

if [ ! TRUE ]; then
   samples=(${samples_orig[@]})

   combined_variants=""
   for sample in "${samples[@]}"
   do
   #	echo "--variant ${result_dir}/${sample}_raw_variants.g.vcf.gz \ "
   	combined_variants+="--variant ${result_dir}/${sample}_raw_variants_AS.g.vcf.gz "
   done
   
   #echo ${!combined_variants}
   
      #-O ${result_dir}/combined_samples_raw_variants_AS.g.vcf.gz \
   gatk_command="gatk CombineGVCFs \
      -R ${ref_fasta} \
      ${combined_variants} \
      -O ${result_dir}/combined_control_samples_raw_variants_AS.g.vcf.gz \
      -G StandardAnnotation \
      -G AS_StandardAnnotation"
   
   echo "gatk_command: ${gatk_command}"
   
   eval "${gatk_command}"
fi   

# CombineGVCFs where each sample GVCFs were generated WITH allele speific annaotation
# This section for control samples: samples_C1_C2 in sample_arrays.sh
# declare -a samples_C1_C2=(
#    C1-gDNA_047_146
#    C1-cfDNA_094_099
#    Urine-Cell-Pellet-C1_023_170
#    C2-gDNA_035_158
#    C2-cfDNA_082_111
#    Urine-Cell-Pellet-C2_011_182
# )

if [ ! TRUE ]; then
   combined_variants=""
   for sample in "${samples_C1_C2[@]}"
   do
   	combined_variants+="--variant ${result_dir}/${sample}_raw_variants_AS.g.vcf.gz "
   done
   
   echo ${combined_variants}

   gatk_command="gatk CombineGVCFs \
      -R ${ref_fasta} \
      ${combined_variants} \
      -O ${result_dir}/combined_control_samples_raw_variants_AS.g.vcf.gz \
      -G StandardAnnotation \
      -G AS_StandardAnnotation"
   
   echo "gatk_command: ${gatk_command}"
   
   eval "${gatk_command}"
fi   

