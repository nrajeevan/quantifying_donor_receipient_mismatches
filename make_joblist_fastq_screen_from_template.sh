#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
result_dir=${work_dir}/results/fastq_screen
out_dir=${work_dir}/FastqScreen

. ${work_dir}/sample_arrays.sh

templateText="module load miniconda/22.11.1; module load BWA/0.7.17-GCCcore-12.2.0; conda activate /gpfs/gibbs/pi/ycga/mane/ycga_bioinfo/envs/fastqScreen; cd ${work_dir}; ${work_dir}/FastqScreen/fastq_screen_xxxxxx.sh";

if [ TRUE ];then

   # 1. All original samples
   samples=(${samples_orig[@]})

   rm -f ${work_dir}/joblist_fastq_screen.txt

   for sample in "${samples[@]}"
   do
      joblist_entry="${templateText/xxxxxx/"${sample}"}";
      echo "${joblist_entry}" >>${work_dir}/joblist_fastq_screen.txt
   done
fi
