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

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/sample_arrays.sh"


templateText="module load FastQC/0.12.1-Java-11; cd ${work_dir}; bash ${work_dir}/FastQC/fastqc_xxxxxx.sh";

if [ TRUE ];then

   # 1. All original samples
   #samples=(${samples_orig[@]})
   samples=(${samples_control[@]})

   rm -f ${work_dir}/joblist_fastqc.txt

   for ((i=0;i< ${#samples[@]} ;i++));
   do
      echo "i: $i"
      sample=${samples[i]}
      echo "${sample}"

      joblist_entry="${templateText/xxxxxx/"${sample}"}";
      echo "${joblist_entry}"
      echo "${work_dir}/joblist_fastqc.txt"
      echo "${joblist_entry}" >>${work_dir}/joblist_fastqc.txt
   done
fi
