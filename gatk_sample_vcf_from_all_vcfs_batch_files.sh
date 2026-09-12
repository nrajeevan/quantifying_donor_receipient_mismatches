#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

scripts_dir=${work_dir}/scripts/GATK

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"

. ${work_dir}/sample_arrays.sh

if [ TRUE ];then

   samples=(${samples_orig[@]})

   for ((i=0;i< ${#samples[@]} ;i++));
   do
      sample=${samples[i]}
      echo "samp1: ${sample}"

      templateFile=${scripts_dir}/gatk_sample_vcf_from_all_vcfs_template.sh
      outfile=${out_dir}/gatk_sample_vcf_from_all_vcfs_${sample}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss/'"${sample}"'/g' $outfile
   done
fi
