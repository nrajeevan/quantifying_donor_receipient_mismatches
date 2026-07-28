#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/fastqc

sample=ssssss
chrom=cccccc

# original sample R1 R2 files
data_dir=${home_dir}/Data/Sample_${sample}

# bbduk trimmed sample R1 R2 files
#data_dir=${work_dir}/results/bbduk

R1_file=${data_dir}/${sample}*R1*
R2_file=${data_dir}/${sample}*R2*

ls ${R1_file}
ls ${R2_file}

if [ ! TRUE ]; then
   chr=${chrom}
   bwa mem -M -t 4 reference_data/hg38/chr${chr} ${R1_file} ${R2_file} >results/bwa/${sample}_chr${chr}.sam
fi

if [ TRUE ]; then
   fastqc -o ${result_dir} -t 4 ${R1_file} ${R2_file}
fi
