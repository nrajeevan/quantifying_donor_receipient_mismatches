#!/bin/bash

# Written by Nallakkandi Rajeevan, Ph.D.
# Contact: n.rajeevan@yale.edu
# Date: July 27, 2026

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
out_dir=${work_dir}/BBDuk

#DIR="${BASH_SOURCE%/*}"
#if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
#. "$DIR/sample_arrays.sh"

. ${work_dir}/sample_arrays.sh

for sample in "${samples_orig[@]}"
do
   templateFile=bbduk_adapter_trim_template.sh
   outfile=${out_dir}/adapter_trim_${sample}.sh
   cp $templateFile $outfile
   sed -i -e 's/ssssss/'"$sample"'/g' $outfile
done
