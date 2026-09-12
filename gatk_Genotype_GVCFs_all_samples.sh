#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

gatk_script_dir=${work_dir}/scripts/GATK

. ${work_dir}/sample_arrays.sh

samples=(${samples_orig[@]})

# CombineGVCFs of all samples
if [ ! TRUE ]; then
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

# genotype GVCFs from all samples
if [ TRUE ]; then

   gatk --java-options "-Xmx32g" GenotypeGVCFs \
   -R ${ref_fasta} \
   -V ${result_dir}/combined_samples_raw_variants_AS.g.vcf.gz \
   -O ${result_dir}/genotypes_from_all_samples_GVCFs_AS.vcf.gz \
   -G StandardAnnotation \
   -G AS_StandardAnnotation
fi

# genotype GVCFs from all control samples
if [  TRUE ]; then

   gatk --java-options "-Xmx32g" GenotypeGVCFs \
   -R ${ref_fasta} \
   -V ${result_dir}/combined_control_samples_raw_variants_AS.g.vcf.gz \
   -O ${result_dir}/genotypes_from_all_control_samples_GVCFs_AS.vcf.gz \
   -G StandardAnnotation \
   -G AS_StandardAnnotation
fi

