#!/bin/bash

home_dir=${HOME}/KidneyTransplant/Mismatch
work_dir=${home_dir}/Analysis
ref_dir=${work_dir}/reference/hg38
ref_fasta=${ref_dir}/hg38.fa
result_dir=${work_dir}/results/gatk
out_dir=${work_dir}/GATK

gatk_script_dir=${work_dir}/scripts/GATK

. ${work_dir}/sample_arrays.sh

samples=(${samples_orig[@]})

if [ ! TRUE ]; then
   combined_variants=""
   for sample in "${samples[@]}"
   do
   #	echo "--variant ${result_dir}/${sample}_raw_variants.g.vcf.gz \ "
   	combined_variants+="--variant ${result_dir}/${sample}_raw_variants.g.vcf.gz "
   done
   
   #echo ${!combined_variants}
   
   gatk_command="gatk CombineGVCFs \
      -R ${ref_fasta} \
      ${combined_variants} \
      -O ${result_dir}/combined_samples_raw_variants.g.vcf.gz"
   
   echo "gatk_command: ${gatk_command}"
   
   eval "${gatk_command}"
fi

# From gatk_Genotype_GVCFs_all_samples.sh
if [ ! TRUE ]; then
   gatk --java-options "-Xmx16g" GenotypeGVCFs \
   -R ${ref_fasta} \
   -V ${result_dir}/combined_samples_raw_variants_AS.g.vcf.gz \
   -O ${result_dir}/genotypes_from_all_samples_GVCFs_AS.vcf.gz \
   -G StandardAnnotation \
   -G AS_StandardAnnotation
fi

# Hard-filter large cohort calset (of the order of 1000s of unrelated samples) using VariantFiltration.
# Note this filter was not used as less number of samples.
if [ ! TRUE ]; then
   in_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS.vcf.gz
   out_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS_excesshet.vcf.gz   

   gatk --java-options "-Xmx16g -Xms3g" VariantFiltration \
    -V ${in_vcf} \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O ${out_vcf}
fi

# Create sites-only VCF using MakeSitesOnlyVcf. Note: Site-level filtering requires 
# only site-level annotation, thid speeding up analysis.
# This retains only first eight columns.
if [ TRUE ]; then
   in_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS.vcf.gz
   out_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS_sitesonly.vcf.gz

   gatk MakeSitesOnlyVcf \
        -I ${in_vcf} \
        -O ${out_vcf}
fi

# Calculate vQSLOD tranches for indels only using VariantRecalibrator.
# Copy Gold_standard_indel_and_snps from cloud:
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
if [ TRUE ]; then
   in_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS_sitesonly.vcf.gz
   out_indels_recal=${result_dir}/genotypes_from_all_samples_GVCFs_AS_indels.recal
   out_indels_tranches=${result_dir}/genotypes_from_all_samples_GVCFs_AS_indels.tranches

   Mills_gold_standard_indels=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   Mills_gold_standard_snps=${ref_dir}/Mills_and_1000G_gold_standard.snps.hg38.vcf.gz
   Axiom_Exome_Plus_genotypes=${ref_dir}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
   Homo_sapiens_assembly38_dbsnp=${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf

   hapmap_3_3_hg38=${ref_dir}/hapmap_3.3.hg38.vcf.gz
   hg38_1000G_omni2_5=${ref_dir}/1000G_omni2.5.hg38.vcf.gz
   hg38_1000G_phase1_snps_high_confidence=${ref_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz


   # INDELs
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
       -V ${in_vcf} \
       --trust-all-polymorphic \
       -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
       -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
       -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
       -mode INDEL \
       --max-gaussians 4 \
       -resource:mills,known=false,training=true,truth=true,prior=12 ${Mills_gold_standard_indels} \
       -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${Axiom_Exome_Plus_genotypes} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${Homo_sapiens_assembly38_dbsnp} \
       -O ${out_indels_recal} \
       --tranches-file ${out_indels_tranches}

   # SNPs
   
   out_snps_recal=${result_dir}/genotypes_from_all_samples_GVCFs_AS_snps.recal
   out_snps_tranches=${result_dir}/genotypes_from_all_samples_GVCFs_AS_snps.tranches
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
       -V ${in_vcf} \
       --trust-all-polymorphic \
       -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 \
       -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
       -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
       -mode SNP \
       --max-gaussians 6 \
       -resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_3_3_hg38} \
       -resource:omni,known=false,training=true,truth=true,prior=12 ${hg38_1000G_omni2_5} \
       -resource:1000G,known=false,training=true,truth=false,prior=10 ${hg38_1000G_phase1_snps_high_confidence} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${Homo_sapiens_assembly38_dbsnp} \
       -O ${out_snps_recal} \
       --tranches-file ${out_snps_tranches}

   # Filter indels on VQSLOD using ApplyVQSR
   in_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS.vcf.gz
   out_indels_recalibrated_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS_indels_recal.vcf.gz
   out_snps_recalibrated_vcf=${result_dir}/genotypes_from_all_samples_GVCFs_AS_snps_recal.vcf.gz
   gatk --java-options "-Xmx24g -Xms24g" \
       ApplyVQSR \
       -V ${in_vcf} \
       --recal-file ${out_indels_recal} \
       --tranches-file ${out_indels_tranches} \
       --truth-sensitivity-filter-level 99.7 \
       --create-output-variant-index true \
       -mode INDEL \
       -O ${out_indels_recalibrated_vcf}

   # Filter SNPs on VQSLOD using ApplyVQSR

   gatk --java-options "-Xmx24g -Xms24g" \
       ApplyVQSR \
       -V ${out_indels_recalibrated_vcf} \
       --recal-file ${out_snps_recal} \
       --tranches-file ${out_snps_tranches} \
       --truth-sensitivity-filter-level 99.7 \
       --create-output-variant-index true \
       -mode SNP \
       -O ${out_snps_recalibrated_vcf}

   #Note: This produces a SNP-filtered callset. Given the indel-filtered callset,
   #      this results in the final filtered callset.

fi

