#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
data_dir=${home_dir}/Data
results_dir=${work_dir}/results/bbduk

sample=ssssss

bbduk=${HOME}/bbmap/bbduk.sh

sample_dir=${data_dir}/Sample_${sample}
R1_file=${sample_dir}/${sample}*R1*.fastq.gz
R2_file=${sample_dir}/${sample}*R2*.fastq.gz
R1_filename=$(basename ${R1_file})
R2_filename=$(basename ${R2_file})

R1_file=${sample_dir}/${R1_filename}
R2_file=${sample_dir}/${R2_filename}
R1_out=${results_dir}/$(basename $R1_file)
R2_out=${results_dir}/$(basename $R2_file)
R1_outm=${results_dir}/Mismatch/$(basename $R1_file)
R2_outm=${results_dir}/Mismatch/$(basename $R2_file)

if [ TRUE ];then
   ${bbduk} \
      in1=${R1_file} \
      in2=${R2_file} \
      out=${R1_out} \
      out2=${R2_out} \
      outm=${R1_outm} \
      outm2=${R2_outm} \
      literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      ktrim=r k=23 mink=11 hdist=1 tpe tbo
fi
