#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis

sample=ssssss

# bbduk trimmed sample R1 R2 files
data_dir=${work_dir}/results/bbduk

out_dir=${work_dir}/results/fastq_screen/tagged_fastq
R1_file=${data_dir}/${sample}*R1*
R2_file=${data_dir}/${sample}*R2*

ls ${R1_file}
ls ${R2_file}

# fastq screen for contamination from other sources (e.g. mouse)
if [ TRUE ];then
   fastq_screen \
   --tag \
   --aligner bwa \
   --threads 8 \
   --conf ${work_dir}/results/fastq_screen/fastq_screen.conf \
   --outdir ${out_dir} \
   ${R1_file}

   fastq_screen \
   --tag \
   --aligner bwa \
   --threads 8 \
   --conf ${work_dir}/results/fastq_screen/fastq_screen.conf \
   --outdir ${out_dir} \
   ${R2_file}

fi
