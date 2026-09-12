#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

. ${work_dir}/sample_arrays.sh

samples=(${samples_orig[@]})

# IndexFeatureFile_GVCFs_all_samples.sh
if [ TRUE ]; then
   for sample in "${samples[@]}"
   do
   	gatk IndexFeatureFile -I ${result_dir}/${sample}_raw_variants_AS.g.vcf.gz
   done
fi   

