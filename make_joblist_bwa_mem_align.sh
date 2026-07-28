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

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"

. ${work_dir}/sample_arrays.sh


templateText="module load BWA/0.7.17-GCCcore-12.2.0; cd ${work_dir}; bash ${work_dir}/BWA/bwa_mem_align_xxxxxx.sh";

if [ TRUE ];then

   # 1. All original samples
   samples=(${samples_orig[@]})

   rm -f ${work_dir}/joblist_bwa_mem_align.txt

   for ((i=0;i< ${#samples[@]} ;i++));
   do
      echo "i: $i"
      sample=${samples[i]}
      echo "${sample}"

      joblist_entry="${templateText/xxxxxx/"${sample}"}";
      echo "${joblist_entry}" >>${work_dir}/joblist_bwa_mem_align.txt
   done
fi
