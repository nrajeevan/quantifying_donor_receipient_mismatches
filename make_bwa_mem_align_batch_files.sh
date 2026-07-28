#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/bwa
out_dir=${work_dir}/BWA

script_dir=${work_dir}/scripts/BWA

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"
. ${work_dir}/sample_arrays.sh

if [ TRUE ];then

   # 1. All original samples
   samples=(${samples_orig[@]})


   echo "${#samples[@]}"
   for ((i=0;i< ${#samples[@]} ;i++));
   do
      echo "i: $i"
      sample=${samples[i]}
      echo "samp1: ${sample1}"

      templateFile=${script_dir}/bwa_mem_align_template.sh
      outfile=${out_dir}/bwa_mem_align_${sample}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss/'"${sample}"'/g' $outfile
   done
fi
