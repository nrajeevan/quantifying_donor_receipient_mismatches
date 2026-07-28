#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

. ${work_dir}/sample_arrays.sh

templateText="module load GATK/4.6.0.0-GCCcore-12.2.0-Java-17; cd ${work_dir}; ${work_dir}/GATK/gatk_variant_call_xxxxxx.sh";

if [ TRUE ];then

   # 1. All original samples
   samples=(${samples_orig[@]})

   rm -f ${work_dir}/joblist_gatk_variant_call.txt

   for sample in "${samples[@]}"
   do
      joblist_entry="${templateText/xxxxxx/"${sample}"}";
      echo "${joblist_entry}" >>${work_dir}/joblist_gatk_variant_call.txt
   done
fi
