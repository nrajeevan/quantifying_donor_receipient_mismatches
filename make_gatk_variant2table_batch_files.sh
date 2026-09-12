#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

scripts_dir=${work_dir}/scripts/GATK

. ${work_dir}/sample_arrays.sh

if [ TRUE ]; then
   samples=(${samples_orig[@]})
   for sample in "${samples[@]}"
   do
      echo "${sample}"
      templateFile=${scripts_dir}/gatk_variant2table_template.sh
      outfile=${out_dir}/gatk_variant2table_${sample}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss/'"${sample}"'/g' $outfile
   done
fi

