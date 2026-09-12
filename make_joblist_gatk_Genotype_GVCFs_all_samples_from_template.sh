#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
results_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

. ${work_dir}/sample_arrays.sh

templateText="module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17; cd ${work_dir}; bash ${work_dir}/GATK/gatk_Genotype_GVCFs_all_samples.sh"

if [ TRUE ];then

   rm -f ${work_dir}/joblist_gatk_Genotype_GVCFs_all_samples.txt
   joblist_entry="${templateText}";
   echo "${joblist_entry}" >>${work_dir}/joblist_gatk_Genotype_GVCFs_all_samples.txt
fi
