#!/bin/bash

sample_id=001

echo "cfDNA: "

echo "AA_AB_AA|AA_BB_AA"
grep -P 'AA_AB_AA|AA_BB_AA' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "AA_AB_AB"
grep -P 'AA_AB_AB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "AA_BB_BB"
grep -P 'AA_BB_BB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "BB_AB_BB|BB_AA_BB"
grep -P 'BB_AB_BB|BB_AA_BB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "BB_AB_AB"
grep -P 'BB_AB_AB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "BB_AA_AA"
grep -P 'BB_AA_AA' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_cfDNA_Rec.txt | wc -l

echo "urDNA: "

echo "AA_AB_AA|AA_BB_AA"
grep -P 'AA_AB_AA|AA_BB_AA' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l

echo "AA_AB_AB"
grep -P 'AA_AB_AB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l

echo "AA_BB_BB"
grep -P 'AA_BB_BB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l

echo "BB_AB_BB|BB_AA_BB"
grep -P 'BB_AB_BB|BB_AA_BB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l

echo "BB_AB_AB"
grep -P 'BB_AB_AB' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l

echo "BB_AA_AA"
grep -P 'BB_AA_AA' Results/sample_${sample_id}_Mismatch_gDNA_Rec_gDNA_Donor_urDNA_Rec.txt | wc -l


