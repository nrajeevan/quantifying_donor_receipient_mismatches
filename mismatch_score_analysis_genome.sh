#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
#home_dir=/Users/n.rajeevan/Madhav/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis

min_depth=10

if [ ! TRUE ];then
   awk -F'\t' '{if($11 != null && $11 != "NA" && ($6 != "INDEL" && $6 != "SNP") && $11 > 20) print $6 "\t" $10 "\t" $11}' CFD-001-cfDNA_092_101_raw_snps.tsv | wc -l
   awk -F'\t' '{if($11 != "NA" && ($6 != "INDEL" && $6 != "SNP") && $11 > 20) print $6 "\t" $10 "\t" $11}' CFD-001-cfDNA_092_101_raw_snps.tsv | wc -l
fi

# For plotting number of variants for different minimum depth
if [ ! TRUE ];then
   for min_depth in {1..50..1}
   do
      #echo "min_depth: $min_depth"
      #1       2       3       4       5       6       7       8       9       10        11        12        13        14
      #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

      if [ ! TRUE ];then
         # Number of variant as a function of minimum depth.
         variant=$(cat CFD-001-cfDNA_092_101_raw_snps.tsv | awk -F'\t' -v depth="${min_depth}" '{if($11 != "NA" && $11 >= depth) print $6 "\t" $10 "\t" $11}' - | wc -l )
         echo "$min_depth	${variant}"
      fi
   done
fi

# All sample IDs

declare -a sample_ids=(
   001
   002
   003
   004
   005
   006
   007
   008
   009
   010
   011
)

if [ TRUE ];then
   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3
fi

find_true_mismatch_detected_in_cfDNA_or_urDNA() {
   local file1=$1
   local file2=$2
   local file3=$3
   local outfile=$4

   echo "${file1}"
   echo "${file2}"
   echo "${file3}"
   echo "${outfile}"

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   if [ ! -f "$file1" ];then
      continue;
   fi

   if [ ! -f "$file2" ];then
      continue;
   fi

   if [ ! -f "$file3" ];then
      continue;
   fi

   if [ -f ${outfile} ];then 
      rm -f ${outfile}; 
   fi

      #awk -F'\t' depth="${min_depth}" -v OFS="\t" -v IFS="," 'BEGIN{First_NR=0; Second_NR=0;}
   if [ TRUE ];then
      awk -F'\t' -v depth="${min_depth}" -v OFS="\t" -v outfile="${outfile}" \
      'BEGIN{First_NR=0; Second_NR=0;}
       NR>1&&NR==FNR{
          if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
   
            key=$1":"$2;
   	       value=$7";"$8";"$10;
   
            gFirst[key]=value;
   
            split($8, ALT, ","); 
            split($10, GT, "[/|]"); 
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) { 
               gFirst_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { 
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gFirst_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gFirst_BB[key]=value;
               }
            }
         }
	 First_NR++;
	 Second_NR++;
         next;
      }
      FNR>1&&NR==(First_NR+FNR+1){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
            gSecond[key]=value;
   
            split($8, ALT, ","); 
            split($10, GT, "[/|]"); 
   
            if($10 == $7"/"$7  || $10 == $7"|"$7) { 
               gSecond_AA[key]=value;
            } 
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { 
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gSecond_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gSecond_BB[key]=value;
               }
            }
         }
         Second_NR++;
         next;
      }
      FNR>1&&NR==(Second_NR+FNR+2){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
            gThird[key]=value;
   
            split($8, ALT, ","); 
            split($10, GT, "[/|]"); 
   
            if($10 == $7"/"$7  || $10 == $7"|"$7) { 
               gThird_AA[key]=value;
            } 
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { 
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gThird_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gThird_BB[key]=value;
               }
            }
         }
      }
      END{
         # file1 to  file2
         if(1) {
            for(key in gFirst_AA) {
               if(key in gSecond_BB) {
                  First_AA_Second_BB[key] = gFirst_AA[key]"_"gSecond_BB[key];

		  if(key in gThird_AA) {
                     First_AA_Second_BB_Third_AA[key] = gFirst_AA[key]"_"gSecond_BB[key]"_"gThird_AA[key];
                  }
                  else if(key in gThird_AB) {
                     First_AA_Second_BB_Third_AB[key] = gFirst_AA[key]"_"gSecond_BB[key]"_"gThird_AB[key];
                  }
                  else if(key in gThird_BB) {
                     First_AA_Second_BB_Third_BB[key] = gFirst_AA[key]"_"gSecond_BB[key]"_"gThird_BB[key];
                  }
                  else {
                     First_AA_Second_BB_Third_None[key] = gFirst_AA[key]"_"gSecond_BB[key]"_None";
                  }
               }
            }

            for(key in gFirst_BB) {
               if(key in gSecond_AA) {
                  First_BB_Second_AA[key] = gFirst_BB[key]"_"gSecond_AA[key];

		  if(key in gThird_AA) {
                     First_BB_Second_AA_Third_AA[key] = gFirst_BB[key]"_"gSecond_AA[key]"_"gThird_AA[key];
                  }
                  else if(key in gThird_AB) {
                     First_BB_Second_AA_Third_AB[key] = gFirst_BB[key]"_"gSecond_AA[key]"_"gThird_AB[key];
                  }
                  else if(key in gThird_BB) {
                     First_BB_Second_AA_Third_BB[key] = gFirst_BB[key]"_"gSecond_AA[key]"_"gThird_BB[key];
                  }
                  else {
                     First_BB_Second_AA_Third_None[key] = gFirst_BB[key]"_"gSecond_AA[key]"_None";
                  }
               }
            }
         }

         # file2 to file1 - removed as same as file1 to  file2
         if(0) {
            for(key in gSecond_AA) {
               if(key in gFirst_BB) {
                  Second_AA_First_BB[key] = gFirst_BB[key]"_"gSecond_AA[key];
               }
            }

            for(key in gSecond_BB) {
               if(key in gFirst_AA) {
                  Second_BB_First_AA[key] = gFirst_AA[key]"_"gSecond_BB[key];
               }
            }
         }

         if(0) {
            for(x in First_AA_Second_BB) print "First_AA_Second_BB", x, First_AA_Second_BB[x] > outfile;
            for(x in First_BB_Second_AA) print "First_BB_Second_AA", x, First_BB_Second_AA[x] > outfile;
         }

	 if(1) {
            for(x in First_AA_Second_BB_Third_AA) print "First_AA_Second_BB_Third_AA", x, First_AA_Second_BB_Third_AA[x] > outfile;
            for(x in First_AA_Second_BB_Third_AB) print "First_AA_Second_BB_Third_AB", x, First_AA_Second_BB_Third_AB[x] > outfile;
            for(x in First_AA_Second_BB_Third_BB) print "First_AA_Second_BB_Third_BB", x, First_AA_Second_BB_Third_BB[x] > outfile;
            for(x in First_AA_Second_BB_Third_None) print "First_AA_Second_BB_Third_None", x, First_AA_Second_BB_Third_None[x] > outfile;

            for(x in First_BB_Second_AA_Third_AA) print "First_BB_Second_AA_Third_AA", x, First_BB_Second_AA_Third_AA[x] > outfile;
            for(x in First_BB_Second_AA_Third_AB) print "First_BB_Second_AA_Third_AB", x, First_BB_Second_AA_Third_AB[x] > outfile;
            for(x in First_BB_Second_AA_Third_BB) print "First_BB_Second_AA_Third_BB", x, First_BB_Second_AA_Third_BB[x] > outfile;
            for(x in First_BB_Second_AA_Third_None) print "First_BB_Second_AA_Third_None", x, First_BB_Second_AA_Third_None[x] > outfile;
         }
      }' "${file1}" "${file2}" "${file3}"
   fi
}

# Run this for generating the "Firest_Second_Third" file
if [ ! TRUE ];then

   suffix=_raw_snps.tsv
   . ${work_dir}/sample_raw_files.sh
   out_dir=${work_dir}/results/gatk/VCF_Tables/MismatchScores

   for sample_id in "${sample_ids[@]}"
   do
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      filePrefix1=gDNA_Donor
      filePrefix2=gDNA_Rec
      filePrefix3=cfDNA_Rec
      #filePrefix3=urDNA_Rec

      echo "out_dir: ${out_dir}"
      outfile="${out_dir}/sample_${sample_id}_Mismatch_Score_${filePrefix1}_${filePrefix2}_${filePrefix3}.txt";

      file1="${samples[$filePrefix1]}${suffix}"
      file2="${samples[$filePrefix2]}${suffix}"
      file3="${samples[$filePrefix3]}${suffix}"

      echo "1. find_mismatch_scores"
      find_true_mismatch_detected_in_cfDNA_or_urDNA ${file1}  ${file2} ${file3} ${outfile}
      echo "2. find_mismatch_scores"
   done
fi

# Run this for generating the output for filling the word document.
if [ TRUE ];then

   . ${work_dir}/sample_raw_files.sh
   out_dir=${work_dir}/results/gatk/VCF_Tables/MismatchScores

   for sample_id in "${sample_ids[@]}"
   do
      echo "Sample ID: ${sample_id}"
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      filePrefix1=gDNA_Donor
      filePrefix2=gDNA_Rec
      filePrefix3=cfDNA_Rec
      #filePrefix3=urDNA_Rec

      outfile="${out_dir}/sample_${sample_id}_Mismatch_Score_${filePrefix1}_${filePrefix2}_${filePrefix3}.txt";
      echo "outfile: ${outfile}"

      n_First_AA_Second_BB=`grep 'First_AA_Second_BB' "${outfile}" | wc -l`
      n_First_AA_Second_BB_Third_AA=`grep 'First_AA_Second_BB_Third_AA' "${outfile}" | wc -l`
      n_First_AA_Second_BB_Third_AB=`grep 'First_AA_Second_BB_Third_AB' "${outfile}" | wc -l`
      n_First_AA_Second_BB_Third_BB=`grep 'First_AA_Second_BB_Third_BB' "${outfile}" | wc -l`
      n_First_AA_Second_BB_Third_None=`grep 'First_AA_Second_BB_Third_None' "${outfile}" | wc -l`


      n_First_BB_Second_AA=`grep 'First_BB_Second_AA' "${outfile}" | wc -l`
      n_First_BB_Second_AA_Third_AA=`grep 'First_BB_Second_AA_Third_AA' "${outfile}" | wc -l`
      n_First_BB_Second_AA_Third_AB=`grep 'First_BB_Second_AA_Third_AB' "${outfile}" | wc -l`
      n_First_BB_Second_AA_Third_BB=`grep 'First_BB_Second_AA_Third_BB' "${outfile}" | wc -l`
      n_First_BB_Second_AA_Third_None=`grep 'First_BB_Second_AA_Third_None' "${outfile}" | wc -l`

      echo "${n_First_AA_Second_BB}	${n_First_BB_Second_AA}"
      echo "n_First_AA_Second_BB:	${n_First_AA_Second_BB_Third_AA}	${n_First_AA_Second_BB_Third_AB}	${n_First_AA_Second_BB_Third_BB}	${n_First_AA_Second_BB_Third_None}	n_First_BB_Second_AA:	${n_First_BB_Second_AA_Third_AA}	${n_First_BB_Second_AA_Third_AB}	${n_First_BB_Second_AA_Third_BB}	${n_First_BB_Second_AA_Third_None}"
   done
fi


find_mismatch_scores() {
   local file1=$1
   local file2=$2
   local outfile=$3

   echo "${file2}"
   echo "${file1}"

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ
      if [ ! -f "$file1" ];then
         continue;
      fi

      if [ ! -f "$file2" ];then
         continue;
      fi

      echo "${outfile}"
      if [ -f ${outfile} ]; then rm -f ${outfile}; fi

      echo "${outfile}"
      if [ TRUE ];then
         awk -F'\t' -v depth="${min_depth}" -v out_file="${outfile}" -v OFS="\t" -v IFS="," \
            'BEGIN{First_NR=0;}
   	 NR>1&&NR==FNR{
            if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
      
               key=$1":"$2;
      	    value=$7";"$8";"$10;
      
               gFirst[key]=value;
      
               split($8, ALT, ","); \
               split($10, GT, "[/|]"); \
      
               if($10 == $7"/"$7 || $10 == $7"|"$7) { 
                  gFirst_AA[key]=value;
               } \
               else {
                  nAlt=0;
                  for(i=1;i<= length(ALT);i++) { \
                     for(j=1;j<=length(GT);j++) {
                        if(GT[j] != "*" && GT[j] == ALT[i]) {
                           nAlt+=1;
                        }
                     }
                  }
      
                  if(nAlt == 1) {
                     gFirst_AB[key]=value;
                  }
      
                  if(nAlt == 2) {
                     gFirst_BB[key]=value;
                  }
               }
            }
   	 First_NR+=1;
            next;
         }
         FNR>1&&NR==(First_NR+FNR+1){
            if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
               key=$1":"$2;
               value=$7";"$8";"$10;
      
               gSecond[key]=value;
      
               split($8, ALT, ","); \
               split($10, GT, "[/|]"); \
      
               if($10 == $7"/"$7  || $10 == $7"|"$7) { 
                  gSecond_AA[key]=value;
               } \
               else {
                  nAlt=0;
                  for(i=1;i<= length(ALT);i++) { \
                     for(j=1;j<=length(GT);j++) {
                        if(GT[j] != "*" && GT[j] == ALT[i]) {
                           nAlt+=1;
                        }
                     }
                  }
      
                  if(nAlt == 1) {
                     gSecond_AB[key]=value;
                  }
      
                  if(nAlt == 2) {
                     gSecond_BB[key]=value;
                  }
               }
            }
      
            if(0) {
               if(key in gFirst_AA) {
                  if(gFirst_AB[key] == value) {
                     print key, gFirst_AB[key], value;
                  }
               }
            }
      
            if(0) {
               if(key in gFirst_BB) {
                  if(gFirst_AB[key] == value) {
                     print key, gFirst_AB[key], value;
                  }
               }
            }
      
            if(0) {
               if(key in gFirst_BB) {
                  if(gFirst_BB[key] == value) {
                     print key, gFirst_BB[key], value;
                  }
               }
            }
         }
         END{
            # file1 to  file2
            if(1) {
               for(key in gFirst_AA) {
                  if(key in gSecond_AA) {
   		     print "AA_AA: 1-2", key, gFirst_AA[key], gSecond_AA[key] >> out_file;
                  }
   	       
                  if(key in gSecond_AB) {
   		     print "AA_AB: 1-2", key, gFirst_AA[key], gSecond_AB[key] >> out_file;
                  }
   	       
   	          if(key in gSecond_BB) {
   		     print "AA_BB: 1-2", key, gFirst_AA[key], gSecond_BB[key] >> out_file;
                  }
      	       }
            }
      
            if(1) {
               for(key in gFirst_BB) {
                  if(key in gSecond_BB) {
   		     print "BB_BB: 1-2", key, gFirst_BB[key], gSecond_BB[key] >> out_file;
                  }
   
                  if(key in gSecond_AB) {
   		     print "BB_AB: 1-2", key, gFirst_BB[key], gSecond_AB[key] >> out_file;
                  }
   
                  if(key in gSecond_AA) {
   		     print "BB_AA: 1-2", key, gFirst_BB[key], gSecond_AA[key] >> out_file;
                  }
               }
            }

            # file2 to file1
            if(1) {
               for(key in gSecond_AA) {
                  if(key in gFirst_AA) {
   		     print "AA_AA: 2-1", key, gFirst_AA[key], gSecond_AA[key] >> out_file;
                  }
   	       
                  if(key in gFirst_AB) {
   		     print "AB_AA: 2-1", key, gFirst_AB[key], gSecond_AA[key] >> out_file;
                  }
   	       
   	          if(key in gFirst_BB) {
   		     print "BB_AA: 2-1", key, gFirst_BB[key], gSecond_AA[key] >> out_file;
                  }
      	       }
            }
      
            if(1) {
               for(key in gSecond_BB) {
                  if(key in gFirst_BB) {
   		     print "BB_BB: 2-1", key, gFirst_BB[key], gSecond_BB[key] >> out_file;
                  }
   
                  if(key in gFirst_AB) {
   		     print "AB_BB: 2-1", key, gFirst_AB[key], gSecond_BB[key] >> out_file;
                  }
   
                  if(key in gFirst_AA) {
   		     print "AA_BB: 2-1", key, gFirst_AA[key], gSecond_BB[key] >> out_file;
                  }
               }
            }

         }' "${file1}" "${file2}"
   fi
}

find_mismatch_scores_old() {
   local file1=$1
   local file2=$2

   echo "${file2}"
   echo "${file1}"

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   . ${work_dir}/sample_raw_files.sh

   out_dir=${work_dir}/results/gatk/VCF_Tables/MismatchScores
   suffix=_raw_snps.tsv

   for sample_id in "${sample_ids[@]}"
   do
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      #file1=gDNA_Rec
      #file2=gDNA_Donor

      file1=gDNA_Donor
      file2=urDNA_Rec

      out_file="${out_dir}/sample_${sample_id}_Mismatch_${file1}_${file2}.txt";

      echo ${file1}
      echo ${file2}
      echo "${samples[$file1]}${suffix}"
      echo "${samples[$file2]}${suffix}"
      echo "${out_file}"
      echo "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"
      echo ""

      if [ TRUE ];then
         continue;
      fi

      if [ ! -f "${samples[$file1]}${suffix}" ];then
         continue;
      fi

      if [ ! -f "${samples[$file2]}${suffix}" ];then
         continue;
      fi

      if [ -f ${out_file} ]; then rm ${out_file}; fi

      awk -F'\t' -v depth="${min_depth}" -v out_file="${out_file}" -v OFS="\t" -v IFS="," \
         'BEGIN{First_NR=0;}
	 NR>1&&NR==FNR{
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
   
            key=$1":"$2;
   	    value=$7";"$8";"$10;
   
            gFirst[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) { 
               gFirst_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gFirst_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gFirst_BB[key]=value;
               }
            }
         }
	 First_NR+=1;
         next;
      }
      FNR>1&&NR==(First_NR+FNR+1){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
            gSecond[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7  || $10 == $7"|"$7) { 
               gSecond_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gSecond_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gSecond_BB[key]=value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_AA) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_BB[key] == value) {
                  print key, gFirst_BB[key], value;
               }
            }
         }
   
      }
      END{
         if(1) {
            for(key in gFirst_AA) {
               if(key in gSecond_AA) {
		  print "AA_AA: ", key, gFirst_AA[key], gSecond_AA[key] >> out_file;
               }
	       
               if(key in gSecond_AB) {
		  print "AA_AB: ", key, gFirst_AA[key], gSecond_AB[key] >> out_file;
               }
	       
	       if(key in gSecond_BB) {
		  print "AA_BB: ", key, gFirst_AA[key], gSecond_BB[key] >> out_file;
               }

   	    }
         }
   
         if(1) {
            for(key in gFirst_BB) {
               if(key in gSecond_BB) {
		  print "BB_BB: ", key, gFirst_BB[key], gSecond_BB[key] >> out_file;
               }

               if(key in gSecond_AB) {
		  print "BB_AB: ", key, gFirst_BB[key], gSecond_AB[key] >> out_file;
               }

               if(key in gSecond_AA) {
		  print "BB_AA: ", key, gFirst_BB[key], gSecond_AA[key] >> out_file;
               }

            }
         }
      }' "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"
   done
}

if [ ! TRUE ];then

   suffix=_raw_snps.tsv
   . ${work_dir}/sample_raw_files.sh
   out_dir=${work_dir}/results/gatk/VCF_Tables/MismatchScores

   for sample_id in "${sample_ids[@]}"
   do
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      filePrefix1=gDNA_Donor
      filePrefix2=gDNA_Rec

      echo "out_dir: ${out_dir}"
      outfile="${out_dir}/sample_${sample_id}_Mismatch_Score_${filePrefix1}_${filePrefix2}.txt";

      file1="${samples[$filePrefix1]}${suffix}"
      file2="${samples[$filePrefix2]}${suffix}"
      #echo "${samples[$filePrefix1]}${suffix}"
      #echo "${samples[$filePrefix2]}${suffix}"
      echo "file1: ${file1}"
      echo "file2: ${file2}"
      echo "${out_file}"
      echo ""

      echo "1. find_mismatch_scores"
      find_mismatch_scores ${file1}  ${file2} ${outfile}
      echo "2. find_mismatch_scores"
   done
fi

# For computing number of AA, AB, BB variants in VCF files (converted to tsv files using GATK VariantsToTable function
if [ ! TRUE ];then

   . ${work_dir}/sample_raw_files.sh

   sample_id=004

   sample_array=sample_${sample_id}[@]
   samples=(${!sample_array})

   for sample in "${samples[@]}"
   do
      echo "sample:  ${sample}"
      if [[ "${sample}" == "none" ]];then
         echo "skipping";
         continue;
      fi
      suffix=_raw_snps.tsv

      if [ ! TRUE ];then
         echo "CFD-001: "
   
         declare -a raw_files=(
            CFD-001_gDNA-Donor_raw_snps.tsv
            CFD-001_gDNA-Rec_raw_snps.tsv
            CFD-001-cfDNA-Rec_raw_snps.tsv
            CFD-001-cfDNA_092_101_raw_snps.tsv
         )
      fi
   
      if [ ! TRUE ];then
         declare -a raw_files=(
            CFD-001_gDNA-Donor_raw_snps.tsv
            CFD-001_gDNA-Rec_raw_snps.tsv
            CFD-001-cfDNA-Rec_raw_snps.tsv
            CFD-001-cfDNA_092_101_raw_snps.tsv
         )
      fi



      #for ((i=0;i< ${#raw_files[@]} ;i++));
      #do
      #file=${raw_files[i]}
      #echo $file

      raw_file=${sample}${suffix}
      cat $raw_file | \
         awk -F'\t' -v depth="${min_depth}" -v OFS="\t" -v IFS="," 'NR>1{ \
	       if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") { \
		  split($8, ALT, ","); \
		  split($10, GT, "[/|]"); \
               
               if($10 == $7"/"$7) { AA+=1;} \

		  nAlt=0;
		  for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
	          }
		  if(nAlt == 1) AB+=1;
		  if(nAlt == 2) BB+=1;

            } \
         } END{
	    print "AA: " AA, "AB: " AB, "BB: " BB, "AA+AB+BB: ", AA+AB+BB;  
	    print AA, AB, BB, AA+AB+BB;} ' -  
   done
fi

# For all patient, find AA, AB, or BB in gDNA-Donor or gDNA-Rec detected in cfDNA or urDNA
if [ ! TRUE ];then

   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   . ${work_dir}/sample_raw_files.sh

   out_dir=Results
   suffix=_raw_snps.tsv

   for sample_id in "${sample_ids[@]}"
   do
      #sample_id=001
   
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})
   
      #for sample in "${samples[@]}"
      #do

      echo "${samples[gDNA_Donor]}"

      file1=gDNA_Rec
      file2=urDNA_Rec

      out_file=${out_dir}/sample_${sample_id}_${file1}_in_${file2};
      echo ${file1}
      echo "${samples[$file1]}${suffix}"
      echo "${samples[$file2]}${suffix}"
      echo "${out_file}"

      echo "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"

      if [ ! TRUE ];then
         continue;
      fi

      awk -F'\t' -v depth="${min_depth}" -v out_file="${out_file}" -v OFS="\t" -v IFS="," \
         'BEGIN{First_NR=0; AA_file=out_file"_AA.txt"; print AA_file; tee_AA = "tee AA_file"; print tee_AA;}
          NR>1&&NR==FNR{ 
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
   
            key=$1":"$2;
   	    value=$7";"$8";"$10;
   
            gFirst[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7) { 
               gFirst_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gFirst_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gFirst_BB[key]=value;
               }
            }
         }
	 First_NR+=1;
         next;
      }
      FNR>1&&NR==(First_NR+FNR+1){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
      
            if(1) {
               if(key in gFirst_AA) {
                  if(gFirst_AA[key] == value) {
                     print key, gFirst_AA[key], value >> out_file"_AA.txt";
                  }
               }
            }
      
            if(1) {
               if(key in gFirst_AB) {
                  if(gFirst_AB[key] == value) {
                     print key, gFirst_AB[key], value >> out_file"_AB.txt";
                  }
               }
            }
      
            if(1) {
               if(key in gFirst_BB) {
                  if(gFirst_BB[key] == value) {
                     print key, gFirst_BB[key], value >> out_file"_BB.txt";
                  }
               }
            }
         }
      }' "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"
   done
fi

# For patient CFD-001 find:  AA gDNA-Rec and (AB or BB) in gDNA-Donor or BB in gDNA-Rec and (AB or AA) in gDNA-Donor
if [ ! TRUE ];then
   declare -a raw_files=(
      CFD-001_gDNA-Donor_raw_snps.tsv
      CFD-001_gDNA-Rec_raw_snps.tsv
      CFD-001-cfDNA-Rec_raw_snps.tsv
      CFD-001-cfDNA_092_101_raw_snps.tsv
   )

   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   awk -F'\t' -v depth="${min_depth}" -v OFS="\t" -v IFS="," \
      'BEGIN{First_NR=0;}
      NR>1&&NR==FNR{ \
      if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {

         key=$1":"$2;
	 value=$7";"$8";"$10;

         gFirst[key]=value;

         split($8, ALT, ","); \
         split($10, GT, "[/|]"); \

         if($10 == $7"/"$7) { 
            gFirst_AA[key]=value;
         } \
         else {
            nAlt=0;
            for(i=1;i<= length(ALT);i++) { \
               for(j=1;j<=length(GT);j++) {
                  if(GT[j] != "*" && GT[j] == ALT[i]) {
                     nAlt+=1;
                  }
               }
            }

            if(nAlt == 1) {
               gFirst_AB[key]=value;
            }

            if(nAlt == 2) {
               gFirst_BB[key]=value;
            }
         }
      }
      First_NR+=1;
      next;
   }
   FNR>1&&NR==(First_NR+FNR+1){
      if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
         key=$1":"$2;
         value=$7";"$8";"$10;

         gSecond[key]=value;

         split($8, ALT, ","); \
         split($10, GT, "[/|]"); \

         if($10 == $7"/"$7) { 
            gSecond_AA[key]=value;
         } \
         else {
            nAlt=0;
            for(i=1;i<= length(ALT);i++) { \
               for(j=1;j<=length(GT);j++) {
                  if(GT[j] != "*" && GT[j] == ALT[i]) {
                     nAlt+=1;
                  }
               }
            }

            if(nAlt == 1) {
               gSecond_AB[key]=value;
            }

            if(nAlt == 2) {
               gSecond_BB[key]=value;
            }
         }
      }

      if(0) {
         if(key in gFirst_AA) {
            if(gFirst_AB[key] == value) {
               print key, gFirst_AB[key], value;
            }
         }
      }

      if(0) {
         if(key in gFirst_BB) {
            if(gFirst_AB[key] == value) {
               print key, gFirst_AB[key], value;
            }
         }
      }

      if(0) {
         if(key in gFirst_BB) {
            if(gFirst_BB[key] == value) {
               print key, gFirst_BB[key], value;
            }
         }
      }

   }
   END{
      if(0) {
         for(key in gFirst_AA) {
            if(key in gSecond_AB) {
               print key, gFirst_AA[key], gSecond_AB[key];
            }
	 }
      }

      if(0) {
         for(key in gFirst_AA) {
            if(key in gSecond_BB) {
               print key, gFirst_AA[key], gSecond_BB[key];
            }
	 }
      }

      if(0) {
         for(key in gFirst_BB) {
            if(key in gSecond_AB) {
               print key, gFirst_BB[key], gSecond_AB[key];
            }
         }
      }

      if(1) {
         for(key in gFirst_BB) {
            if(key in gSecond_AA) {
               print key, gFirst_BB[key], gSecond_AA[key];
            }
         }
      }

   }' ${raw_files[gDNA_Rec]} ${raw_files[gDNA_Donor]}
fi

# For All patients find:  AA gDNA-Rec and (AB or BB) in gDNA-Donor or BB in gDNA-Rec and (AB or AA) in gDNA-Donor
# Note: From the out_file use e.g.: grep 'AA_AB:' Results/sample_001_Mismatch_gDNA_Rec_gDNA_Donor.txt | wc -l
if [ ! TRUE ];then

   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   . ${work_dir}/sample_raw_files.sh

   out_dir=Results
   suffix=_raw_snps.tsv

   for sample_id in "${sample_ids[@]}"
   do
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      #file1=gDNA_Rec
      #file2=gDNA_Donor

      file1=gDNA_Donor
      file2=urDNA_Rec

      out_file="${out_dir}/sample_${sample_id}_Mismatch_${file1}_${file2}.txt";

      echo ${file1}
      echo ${file2}
      echo "${samples[$file1]}${suffix}"
      echo "${samples[$file2]}${suffix}"
      echo "${out_file}"
      echo "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"
      echo ""

      if [ ! TRUE ];then
         continue;
      fi

      if [ ! -f "${samples[$file1]}${suffix}" ];then
         continue;
      fi

      if [ ! -f "${samples[$file2]}${suffix}" ];then
         continue;
      fi

      if [ -f ${out_file} ]; then rm ${out_file}; fi

      awk -F'\t' -v depth="${min_depth}" -v out_file="${out_file}" -v OFS="\t" -v IFS="," \
         'BEGIN{First_NR=0;}
	 NR>1&&NR==FNR{
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
   
            key=$1":"$2;
   	    value=$7";"$8";"$10;
   
            gFirst[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) { 
               gFirst_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gFirst_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gFirst_BB[key]=value;
               }
            }
         }
	 First_NR+=1;
         next;
      }
      FNR>1&&NR==(First_NR+FNR+1){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
            gSecond[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7  || $10 == $7"|"$7) { 
               gSecond_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gSecond_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gSecond_BB[key]=value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_AA) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_BB[key] == value) {
                  print key, gFirst_BB[key], value;
               }
            }
         }
   
      }
      END{
         if(1) {
            for(key in gFirst_AA) {
               if(key in gSecond_AA) {
		  print "AA_AA: ", key, gFirst_AA[key], gSecond_AA[key] >> out_file;
               }
	       
               if(key in gSecond_AB) {
		  print "AA_AB: ", key, gFirst_AA[key], gSecond_AB[key] >> out_file;
               }
	       
	       if(key in gSecond_BB) {
		  print "AA_BB: ", key, gFirst_AA[key], gSecond_BB[key] >> out_file;
               }

   	    }
         }
   
         if(1) {
            for(key in gFirst_BB) {
               if(key in gSecond_BB) {
		  print "BB_BB: ", key, gFirst_BB[key], gSecond_BB[key] >> out_file;
               }

               if(key in gSecond_AB) {
		  print "BB_AB: ", key, gFirst_BB[key], gSecond_AB[key] >> out_file;
               }

               if(key in gSecond_AA) {
		  print "BB_AA: ", key, gFirst_BB[key], gSecond_AA[key] >> out_file;
               }

            }
         }
      }' "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}"
   done
fi


# For All patients find:  AA gDNA-Rec and (AB or BB) in gDNA-Donor or BB in gDNA-Rec and (AB or AA) in gDNA-Donor
# Note: From the out_file use e.g.: grep 'AA_AB:' Results/sample_001_Mismatch_gDNA_Rec_gDNA_Donor.txt | wc -l
if [ ! TRUE ];then

   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   . ${work_dir}/sample_raw_files.sh

   out_dir=Results
   suffix=_raw_snps.tsv

   for sample_id in "${sample_ids[@]}"
   do
      sample_array=sample_${sample_id}[@]
      samples=(${!sample_array})

      file1=gDNA_Rec
      file2=gDNA_Donor
      file3=cfDNA_Rec

      out_file="${out_dir}/sample_${sample_id}_Mismatch_${file1}_${file2}_${file3}.txt";

      echo ${file1}
      echo ${file2}
      echo ${file3}
      echo "${samples[$file1]}${suffix}"
      echo "${samples[$file2]}${suffix}"
      echo "${samples[$file3]}${suffix}"
      echo "${out_file}"
      echo "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}" "${samples[$file3]}${suffix}"
      echo ""

      if [ ! TRUE ];then
         continue;
      fi

      if [ ! -f "${samples[$file1]}${suffix}" ]; then
         continue;
      fi

      if [ ! -f "${samples[$file2]}${suffix}" ]; then
         continue;
      fi

      if [ ! -f "${samples[$file3]}${suffix}" ]; then
         continue;
      fi

      if [ -f ${out_file} ]; then rm ${out_file}; fi

      awk -F'\t' -v depth="${min_depth}" -v out_file="${out_file}" -v OFS="\t" -v IFS="," \
      'BEGIN{First_NR=0;Second_NR=0;} NR>1&&NR==FNR{ \
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
   
            key=$1":"$2;
   	    value=$7";"$8";"$10;
   
            gFirst[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) { 
               gFirst_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gFirst_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gFirst_BB[key]=value;
               }
            }
         }
	 First_NR+=1;
	 Second_NR+=1;
         next;
      }
      FNR>1&&NR==(First_NR+FNR+1){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
            gSecond[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) {
               gSecond_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gSecond_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gSecond_BB[key]=value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_AA) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_AB[key] == value) {
                  print key, gFirst_AB[key], value;
               }
            }
         }
   
         if(0) {
            if(key in gFirst_BB) {
               if(gFirst_BB[key] == value) {
                  print key, gFirst_BB[key], value;
               }
            }
         }
         Second_NR+=1;
         next;
      }
      FNR>1&&NR==(Second_NR+FNR+2){
         if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
            key=$1":"$2;
            value=$7";"$8";"$10;
   
	    if(key == "chr3:39143468") {
               print NR, FNR
               print key, value;
            }

            gThird[key]=value;
   
            split($8, ALT, ","); \
            split($10, GT, "[/|]"); \
   
            if($10 == $7"/"$7 || $10 == $7"|"$7) {
               gThird_AA[key]=value;
            } \
            else {
               nAlt=0;
               for(i=1;i<= length(ALT);i++) { \
                  for(j=1;j<=length(GT);j++) {
                     if(GT[j] != "*" && GT[j] == ALT[i]) {
                        nAlt+=1;
                     }
                  }
               }
   
               if(nAlt == 1) {
                  gThird_AB[key]=value;
               }
   
               if(nAlt == 2) {
                  gThird_BB[key]=value;
               }
            }
         }
   
         next;
      }
      END{

         print length(gFirst_AA);
         print length(gFirst_AB);
         print length(gFirst_BB);
         print length(gSecond_AA);
         print length(gSecond_AB);
         print length(gSecond_BB);
         print length(gThird_AA);
         print length(gThird_AB);
         print length(gThird_BB);
         if(1) {
            for(key in gFirst_AA) {
               if(key in gSecond_AB) {
                  if(key in gThird_AA) {
                     if(gFirst_AA[key] == gThird_AA[key]) {
			print "AA_AB_AA: ", key, gFirst_AA[key], gSecond_AB[key], gThird_AA[key] >> out_file;
                     }
                  }
                  if(key in gThird_AB) {
                     if(gSecond_AB[key] == gThird_AB[key]) {
			print "AA_AB_AB: ", key, gFirst_AA[key], gSecond_AB[key], gThird_AB[key] >> out_file;
                     }
                  }
               }

               if(key in gSecond_BB) {
                  if(key in gThird_AA) {
                     if(gFirst_AA[key] == gThird_AA[key]) {
			print "AA_BB_AA: ", key, gFirst_AA[key], gSecond_BB[key], gThird_AA[key] >> out_file;
                     }
                  }

		  if(key == "chr3:39143468") {
                     print key, gFirst_AA[key], gSecond_BB[key], gThird_BB[key];
		  }

                  if(key in gThird_BB) {
                     if(gSecond_BB[key] == gThird_BB[key]) {
			print "AA_BB_BB: ", key, gFirst_AA[key], gSecond_BB[key], gThird_BB[key] >> out_file;
                     }
                  }
               }
   	    }

            for(key in gFirst_BB) {
               if(key in gSecond_AB) {
                  if(key in gThird_BB) {
                     if(gFirst_BB[key] == gThird_BB[key]) {
			print "BB_AB_BB: ", key, gFirst_BB[key], gSecond_AB[key], gThird_BB[key] >> out_file;
                     }
                  }
                  if(key in gThird_AB) {
                     if(gSecond_AB[key] == gThird_AB[key]) {
			print "BB_AB_AB: ", key, gFirst_BB[key], gSecond_AB[key], gThird_AB[key] >> out_file;
                     }
                  }
               }

               if(key in gSecond_AA) {
                  if(key in gThird_BB) {
                     if(gFirst_BB[key] == gThird_BB[key]) {
			print "BB_AA_BB: ", key, gFirst_BB[key], gSecond_AA[key], gThird_BB[key] >> out_file;
                     }
                  }
                  if(key in gThird_AA) {
                     if(gSecond_AA[key] == gThird_AA[key]) {
			print "BB_AA_AA: ", key, gFirst_BB[key], gSecond_AA[key], gThird_AA[key] >> out_file;
                     }
                  }
               }
   	    }
         }
      }' "${samples[$file1]}${suffix}" "${samples[$file2]}${suffix}" "${samples[$file3]}${suffix}"
   done
fi



# To find AA, AB, or BB in gDNA-Rec and detected in cfDNA or urDNA.
if [ ! TRUE ];then
   declare -a raw_files=(
      CFD-001_gDNA-Donor_raw_snps.tsv
      CFD-001_gDNA-Rec_raw_snps.tsv
      CFD-001-cfDNA-Rec_raw_snps.tsv
      CFD-001-cfDNA_092_101_raw_snps.tsv
   )

   gDNA_Donor=0
   gDNA_Rec=1
   cfDNA_Rec=2
   urDNA_Rec=3

   #echo "${raw_files[gDNA_Donor]}, ${raw_files[urDNA_Rec]}";
   #1       2       3       4       5       6       7       8       9       10        11        12        13        14
   #CHROM   POS     ID      QUAL    FILTER  TYPE    REF     ALT     AF      GT        DP        PL        GP        GQ

   awk -F'\t' -v depth="${min_depth}" -v OFS="\t" -v IFS="," \
      'BEGIN{First_NR=0;} NR>1&&NR==FNR{ \
      if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {

         key=$1":"$2;
	 value=$7";"$8";"$10;

         gFirst[key]=value;

         split($8, ALT, ","); \
         split($10, GT, "[/|]"); \

	 if(key == "chr3:36546049") {
            print "gDNA_Rec: ", key, value;
         }

         if($10 == $7"/"$7) { 
            gFirst_AA[key]=value;
         } \
         else {
            nAlt=0;
            for(i=1;i<= length(ALT);i++) { \
               for(j=1;j<=length(GT);j++) {
                  if(GT[j] != "*" && GT[j] == ALT[i]) {
                     nAlt+=1;
                  }
               }
            }

            if(nAlt == 1) {
               gFirst_AB[key]=value;
            }

            if(nAlt == 2) {
               gFirst_BB[key]=value;
            }
         }
      }
      First_NR+=1;
      next;
   }
   FNR>1&&NR==(First_NR+FNR+1){
      if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
         key=$1":"$2;
         value=$7";"$8";"$10;

         gSecond[key]=value;

         split($8, ALT, ","); \
         split($10, GT, "[/|]"); \

	 if(key == "chr3:36546049") {
            print "gDNA_Donor: ", key, value;
         }

         if($10 == $7"/"$7) { 
            gSecond_AA[key]=value;
         } \
         else {
            nAlt=0;
            for(i=1;i<= length(ALT);i++) { \
               for(j=1;j<=length(GT);j++) {
                  if(GT[j] != "*" && GT[j] == ALT[i]) {
                     nAlt+=1;
                  }
               }
            }

            if(nAlt == 1) {
               gSecond_AB[key]=value;
            }

            if(nAlt == 2) {
               gSecond_BB[key]=value;
            }
         }
      }

      next;
   }

   FNR>1{
      if($11 != null && $11 != "NA" && $11 >= depth && $10 != "./." && $10 != "*/*" && $10 != "*|*") {
         key=$1":"$2;
         value=$7";"$8";"$10;

         gThird[key]=value;

         split($8, ALT, ","); \
         split($10, GT, "[/|]"); \

	 if(key == "chr3:36546049") {
            print "cfDNA_Rec: ", key, value;
         }

         if($10 == $7"/"$7) { 
            gThird_AA[key]=value;
         } \
         else {
            nAlt=0;
            for(i=1;i<= length(ALT);i++) { \
               for(j=1;j<=length(GT);j++) {
                  if(GT[j] != "*" && GT[j] == ALT[i]) {
                     nAlt+=1;
                  }
               }
            }

            if(nAlt == 1) {
               gThird_AB[key]=value;
            }

            if(nAlt == 2) {
               gThird_BB[key]=value;
            }
         }
      }

      next;
   }
   END{
      if(0) {
         for(key in gFirst_AA) {
            if(key in gSecond_AB) {
               if(key in gThird_AB) {
                  if(gSecond_AB[key] == gThird_AB[key]) {
                     print key, gFirst_AA[key], gSecond_AB[key], gThird_AB[key];
                  }
               }
            }
	 }
      }

      if(0) {
         for(key in gFirst_AA) {
            if(key in gSecond_AB) {
               if(key in gThird_AA) {
                  if(gFirst_AA[key] == gThird_AA[key]) {
                     print key, gFirst_AA[key], gSecond_AB[key], gThird_AA[key];
                  }
               }
            }
         }
      }

      if(0) {
         for(key in gFirst_AA) {
            if(key in gSecond_BB) {
               if(key in gThird_BB) {
                  if(gSecond_BB[key] == gThird_BB[key]) {
                     print key, gFirst_AA[key], gSecond_BB[key], gThird_BB[key];
                  }
               }
            }
	 }
      }

      if(0) {
         for(key in gFirst_BB) {
            if(key in gSecond_AB) {
               if(key in gThird_AB) {
                  if(gSecond_AB[key] == gThird_AB[key]) {
                     print key, gFirst_BB[key], gSecond_AB[key], gThird_AB[key];
                  }
               }
            }
         }
      }

      if(0) {
         for(key in gFirst_BB) {
            if(key in gSecond_AB) {
               if(key in gThird_BB) {
                  if(gFirst_BB[key] == gThird_BB[key]) {
                     print key, gFirst_BB[key], gSecond_AB[key], gThird_BB[key];
                  }
               }
            }
         }
      }

      if(0) {
         for(key in gFirst_BB) {
            if(key in gSecond_AA) {
               if(key in gThird_AA) {
                  if(gSecond_AA[key] == gThird_AA[key]) {
                     print key, gFirst_BB[key], gSecond_AA[key], gThird_AA[key];
                  }
               }
            }
         }
      }

   }' ${raw_files[gDNA_Rec]} ${raw_files[gDNA_Donor]} ${raw_files[urDNA_Rec]}
fi

