#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

scripts_dir=${work_dir}/scripts/GATK

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"

. ${work_dir}/sample_arrays.sh

if [ TRUE ];then

   # 1. Variants in donor gDNA detected in recipient cfDNA 
   if [ ! TRUE ];then
      samples=(${samples_cfDNA_rec_gDNA_donor[@]})
      input1_ext=raw
      input2_ext=raw
      output_ext=comm
      sel_var_method=concordance
   fi

   # 1. Variants in donor gDNA detected in recipient cfDNA
   if [ TRUE ];then
      samples=(${samples_urDNA_rec_gDNA_donor[@]})
      input1_ext=raw
      input2_ext=raw
      output_ext=comm
      sel_var_method=concordance
   fi

   # 2. Variants in recipient gDNA detected in recipient cfDNA
   if [ ! TRUE ];then
      samples=(${samples_cfDNA_rec_gDNA_rec[@]})
      input1_ext=raw
      input2_ext=raw
      output_ext=comm
      sel_var_method=concordance
   fi

   # 3. Variants common in recipient and donor gDNA that are detected in recipient cfDNA
   if [ ! TRUE ];then
      samples=(${samples_cfDNA_rec_gDNA_rec_and_donor[@]})
      input1_ext=comm
      input2_ext=comm
      output_ext=comm
      sel_var_method=concordance
   fi

   # 4. Variants in recipient gDNA and not in donor gDNA
   if [ ! TRUE ];then
      samples=(${samples_gDNA_rec_not_gDNA_donor[@]})
      input1_ext=raw
      input2_ext=raw
      output_ext=diff
      sel_var_method=discordance
   fi

   # 5. Variants in donor gDNA and not in recipient gDNA
   if [ ! TRUE ];then
      samples=(${samples_gDNA_donor_not_gDNA_rec[@]})
      input1_ext=raw
      input2_ext=raw
      output_ext=diff
      sel_var_method=discordance
   fi

   # 6. Variants in recipient gDNA and not in donor gDNA that are detected in recipient cfDNA
   if [ ! TRUE ];then
      samples=(${samples_cfDNA_rec_gDNA_rec_not_donor[@]})
      input1_ext=raw
      input2_ext=diff
      output_ext=comm
      sel_var_method=concordance
   fi

   # 7. Variants in donor gDNA and not in recipient gDNA that are detected in recipient cfDNA
   if [ ! TRUE ];then
      samples=(${samples_cfDNA_rec_gDNA_donor_not_rec[@]})
      input1_ext=raw
      input2_ext=diff
      output_ext=comm
      sel_var_method=concordance
   fi

   echo "${#samples[@]}"
   echo "${#samples_rec_rec_cfDNA_gDNA[@]}"
   for ((i=0;i< ${#samples[@]} ;i+=3));
   do
      echo "i: $i"
      sample1=${samples[i]}
      sample2=${samples[i+1]}
      sample12=${samples[i+2]}

      echo "samp1: ${sample1}"
      echo "samp2: ${sample2}"
      echo "samp12: ${sample12}"

      templateFile=${scripts_dir}/gatk_compare_vcfs_template.sh
      outfile=${out_dir}/gatk_compare_vcfs_${sample12}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss12/'"${sample12}"'/g' $outfile
      sed -i -e 's/ssssss1/'"${sample1}"'/g' $outfile
      sed -i -e 's/ssssss2/'"${sample2}"'/g' $outfile

      sed -i -e 's/input1ext/'"${input1_ext}"'/g' $outfile
      sed -i -e 's/input2ext/'"${input2_ext}"'/g' $outfile
      sed -i -e 's/outputext/'"${output_ext}"'/g' $outfile
      sed -i -e 's/selectvariant_method/'"${sel_var_method}"'/g' $outfile
   done
fi


DO_SAMPLES_REC_CFDNA_GDNA=TRUE
if [ ! DO_SAMPLES_REC_CFDNA_GDNA ];then
   # .,123s/samples_rec_cfDNA_gDNA/samples_rec_cfDNA_gDNA/g
   echo "${#samples_rec_cfDNA_gDNA[@]}"
   for ((i=0;i< ${#samples_rec_cfDNA_gDNA[@]} ;i+=2));
   do
      echo "i: $i"
      echo "samp1: ${samples_rec_cfDNA_gDNA[i]} samp2: ${samples_rec_cfDNA_gDNA[i+1]}"
      sample1=${samples_rec_cfDNA_gDNA[i]}
      sample2=${samples_rec_cfDNA_gDNA[i+1]}
   
      sample=${sample1%-*}
      sample12=${sample1}_${sample2}
      #echo ${sample}
      templateFile=${scripts_dir}/gatk_compare_vcfs_template.sh
      outfile=${out_dir}/gatk_compare_vcfs_${sample12}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss1/'"${sample1}"'/g' $outfile
      sed -i -e 's/ssssss2/'"${sample2}"'/g' $outfile
   done
fi

DO_SAMPLES_DONOR_CFDNA_GDNA=TRUE
if [ ! DO_SAMPLES_DONOR_CFDNA_GDNA ];then
   echo "${#samples_donor_cfDNA_gDNA[@]}"
   for ((i=0;i< ${#samples_donor_cfDNA_gDNA[@]} ;i+=2));
   do
      echo "i: $i"
      echo "samp1: ${samples_donor_cfDNA_gDNA[i]} samp2: ${samples_donor_cfDNA_gDNA[i+1]}"
      sample1=${samples_donor_cfDNA_gDNA[i]}
      sample2=${samples_donor_cfDNA_gDNA[i+1]}
   
      sample=${sample1%-*}
      sample12=${sample1}_${sample2}
      #echo ${sample}
      templateFile=${scripts_dir}/gatk_compare_vcfs_template.sh
      outfile=${out_dir}/gatk_compare_vcfs_${sample12}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss1/'"${sample1}"'/g' $outfile
      sed -i -e 's/ssssss2/'"${sample2}"'/g' $outfile
   done
fi

# Testing

samples=(${samples_rec_rec_cfDNA_gDNA[@]})
if [ ! TRUE ];then
   echo "${#samples[@]}"
   echo "${#samples_rec_rec_cfDNA_gDNA[@]}"
   for ((i=0;i< ${#samples[@]} ;i+=3));
   do
      echo "i: $i"
      echo "samp1: ${samples[i]}"
      echo "samp2: ${samples[i+1]}"
      echo "samp12: ${samples[i+2]}"
   done
fi

# 1. Variants in donor gDNA detected in recipient cfDNA
DO_SAMPLES_DONOR_REC_CFDNA_GDNA=TRUE
if [ ! DO_SAMPLES_DONOR_REC_CFDNA_GDNA ];then
   echo "${#samples_rec_rec_cfDNA_gDNA[@]}"
   for ((i=0;i< ${#samples_rec_rec_cfDNA_gDNA[@]} ;i+=3));
   do
      echo "i: $i"
      echo "samp1: ${samples_rec_rec_cfDNA_gDNA[i]}"
      echo "samp2: ${samples_rec_rec_cfDNA_gDNA[i+1]}"
      echo "samp12: ${samples_rec_rec_cfDNA_gDNA[i+2]}"
      sample1=${samples_rec_rec_cfDNA_gDNA[i]}
      sample2=${samples_rec_rec_cfDNA_gDNA[i+1]}
      sample12=${samples_rec_rec_cfDNA_gDNA[i+2]}
   
      templateFile=${scripts_dir}/gatk_compare_vcfs_template.sh
      outfile=${out_dir}/gatk_compare_vcfs_${sample12}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss12/'"${sample12}"'/g' $outfile
      sed -i -e 's/ssssss1/'"${sample1}"'/g' $outfile
      sed -i -e 's/ssssss2/'"${sample2}"'/g' $outfile
   done
fi

if [ ! TRUE ];then
   for ((i=0;i< ${#samples[@]} ;i+=2));
   do
      echo "samp1: ${samples[i]} samp2: ${samples[i+1]}"
      sample1=${samples[i]}
      sample2=${samples[i+1]}
   
      sample=${sample1%-*}
      sample12=${sample1}_${sample2}
      echo ${sample}
      templateFile=${scripts_dir}/gatk_compare_vcfs_template.sh
      outfile=${out_dir}/gatk_compare_vcfs_${sample12}.sh
      cp $templateFile $outfile
      sed -i -e 's/ssssss1/'"${sample1}"'/g' $outfile
      sed -i -e 's/ssssss2/'"${sample2}"'/g' $outfile
   done
fi

