#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/fastqc
out_dir=${work_dir}/FastQC

script_dir=${work_dir}/scripts/FastQC

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/sample_arrays.sh"

echo "${DIR}"
echo "${BASH_SOURCE}"
echo "${BASH_SOURCE%/*}"

if [ TRUE ];then

   # 1. All original samples
   #samples=(${samples_orig[@]})
   samples=(${samples_control[@]})


   echo "${#samples[@]}"
   for ((i=0;i< ${#samples[@]} ;i++));
   do
      echo "i: $i"
      sample=${samples[i]}
      echo "samp1: ${sample}"
   if [ TRUE ];then
      templateFile=${script_dir}/fastqc_template.sh
      outfile=${out_dir}/fastqc_${sample}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss/'"${sample}"'/g' $outfile
   fi
   done
fi
