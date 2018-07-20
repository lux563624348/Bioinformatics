#!/bin/bash
#set -e	#### Terminating the script if any command exited with a nonzero exit status
set -u	#### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.
########################################################################
## 8/08/2017
## Last Modified by Xiang Li, basing on Zhouhao's scripts.
## lux@gwu.edu
## Peng's Lab
## Ver.beta
########################################################################

########################################################################
  #  means NEED MODIFICATION. "VERY IMPORTANT INFORMATION"
  ## means Title_level 1.
 ### means Title_level 2.
#### means comment.
########################################################################
echo "Start Date: `date` "
Start_Date=`date`

source ./functions.sh
echo "Import functions.sh" 

#			NEED MODIFICATION FOR DIFFERENT PROJECT

########################################################################
## GLOBAL VARIABLES
########################################################################
__RAW_DATA_PATH_DIR="/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620"
__EXE_PATH="/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620/Overlap_WT_Treg_foxp3_binding"

########################################################################
#echo "__INPUT_SAMPLE_DIR_List=(1_input_Bruce4 2_H3K4me3_Bruce4 3_H3K9me3_Bruce4 4_input_WT 5_H3K4me3_WT 6_input_mir34bc_KO 7_H3K4me3_mir34bc_KO 8_H3K9me3_WT_Dec 8_H3K9me3_WT_April 9_H3K9me3_mir34bc_KO_Dec 9_H3K9me3_mir34bc_KO_April 10_New_H3K9me3_WT 11_New_H3K9me3_mir34bc_KO)"
__INPUT_SAMPLE_DIR_List=(
Sample_foxp3_peaks/GSE40684_foxp3_peaks.bed
Sample_Treg_TCF1_20160827000/SICER_results/Sample_Treg_TCF1_20160827000-W200-G400-FDR0.0001-islandfiltered.bed
)


DIR_CHECK_CREATE ${__EXE_PATH}

####   Overlap for two input.
Title="WT_Treg_Tcf1_vs_Foxp3"
echo "$Title"
echo ""
echo ""
i=1
j=0
Input_A1=${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[$i]}
Input_B1=${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[$j]}
wc -l ${Input_A1}
wc -l ${Input_B1}
Output_Path="${__EXE_PATH}/Overlap_${Title}_ref_${__INPUT_SAMPLE_DIR_List[$i]: 0 : 18 }"
Output_Path2="${__EXE_PATH}/Specific_${Title}_ref_${__INPUT_SAMPLE_DIR_List[$i]: 0 : 18 }"

echo "bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.bed"
bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.bed
echo "bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.bed"
bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.bed

wc -l ${Output_Path}.bed
wc -l ${Output_Path2}.bed
########################################################################

Title="Foxp3_vs_WT_Treg_Tcf1"
echo "$Title"
echo ""
echo ""
i=0
j=1
Input_A1=${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[$i]}
Input_B1=${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[$j]}
wc -l ${Input_A1}
wc -l ${Input_B1}
Output_Path="${__EXE_PATH}/Overlap_${Title}_ref_${__INPUT_SAMPLE_DIR_List[$i]: 0 : 18 }"
Output_Path2="${__EXE_PATH}/Specific_${Title}_ref_${__INPUT_SAMPLE_DIR_List[$i]: 0 : 18 }"

echo "bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.bed"
bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.bed
echo "bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.bed"
bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.bed

wc -l ${Output_Path}.bed
wc -l ${Output_Path2}.bed

EMAIL_ME
