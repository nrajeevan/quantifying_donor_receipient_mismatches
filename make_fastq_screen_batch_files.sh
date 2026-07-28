#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
result_dir=${work_dir}/results/fastq_screen
out_dir=${work_dir}/FastqScreen

scripts_dir=${work_dir}/scripts/FastqScreen

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

      templateFile=${scripts_dir}/fastq_screen_template.sh
      outfile=${out_dir}/fastq_screen_${sample}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss/'"${sample}"'/g' $outfile
   done
fi
