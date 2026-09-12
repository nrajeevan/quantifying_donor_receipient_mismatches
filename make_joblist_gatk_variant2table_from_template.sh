#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"

. ${work_dir}/sample_arrays.sh

rm -f ${work_dir}/joblist_gatk_variant2table.txt

templateText="module load GATK/4.6.0.0-GCCcore-12.2.0-Java-17; cd ${work_dir}; bash ${work_dir}/GATK/gatk_variant2table_xxxxxx.sh";

samples=(${samples_orig[@]})
#samples=(${samples_cfDNA_rec_gDNA_donor_out[@]})
#samples=(${samples_urDNA_rec_gDNA_donor_out[@]})

for sample in "${samples[@]}"
do
   joblist_entry="${templateText/xxxxxx/"${sample}"}";
   echo "${joblist_entry}" >>${work_dir}/joblist_gatk_variant2table.txt
done
