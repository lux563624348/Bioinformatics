#!/bin/bash
#set -e	#### Terminating the script if any command exited with a nonzero exit status
set -u	#### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.

## 
#	nohup bash run_main.sh 1>default.log 2>default_err.log &
########################################################################
## 12/02/2017
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
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
__RAW_DATA_PATH_DIR=/mnt/c/Users/Xiang/Dropbox/AResearch/version.beta/current/BAMPE_PEAKS/TLE3/TEST
#### Execution or Output directory
__EXE_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620/Peaks/data/tcf1_foxp3
########################################################################
#echo "__INPUT_SAMPLE_DIR_List=(1_input_Bruce4 2_H3K4me3_Bruce4 3_H3K9me3_Bruce4 4_input_WT 5_H3K4me3_WT 6_input_mir34bc_KO 7_H3K4me3_mir34bc_KO 8_H3K9me3_WT_Dec 8_H3K9me3_WT_April 9_H3K9me3_mir34bc_KO_Dec 9_H3K9me3_mir34bc_KO_April 10_New_H3K9me3_WT 11_New_H3K9me3_mir34bc_KO)"


__INPUT_SAMPLE_DIR_List=(
Naive_CD4_Rest      #0
Naive_CD4_Active_40min #1
Naive_CD4_Active_150min #2
Naive_CD4_Active_15h
TCM_CD4_Rest   #4
TCM_CD4_Active_40min
TCM_CD4_Active_150min
TCM_CD4_Active_15h
TEM_CD4_Rest   #8
TEM_CD4_Active_40min
TEM_CD4_Active_150min
TEM_CD4_Active_15h #11
)
__INPUT_SAMPLE_DIR_List=(
Common_peaks_Treg_Foxp3
Solo_peaks_Treg_Foxp3
)


echo "INPUT_SAMPLE_DIR_List= (${__INPUT_SAMPLE_DIR_List[*]})"


#__OUT_SAMPLE_NAME_List=(Unstim_Ctrl2_TCM Unstim_TKO1_TCM_bp Unstim_TKO2_TCM)
SAMPLE_NUM=${#__INPUT_SAMPLE_DIR_List[*]}
########################################################################

########################################################################
##	MAIN BODY
########################################################################

main() {
#### Saving DIR Check and Create
DIR_CHECK_CREATE $__RAW_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH

#### Email Alert
FUNC_CHOOSE_EMAIL_ALERT
Alert_email=$?

echo ""
echo "__FASTQ_DIR_R1 __FASTQ_DIR_R2 are the READS_FULL_DIR FOR ANALYSIS"

for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	#PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	RUN_Venn_Diagram ${__RAW_DATA_PATH_DIR} "bed"
	break
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]}
	
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]}
	
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
done


## SECOND LOOP




	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}"
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

