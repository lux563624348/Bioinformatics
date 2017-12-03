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
__RAW_DATA_PATH_DIR=$(pwd)
__RAW_DATA_PATH_DIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620
#### Execution or Output directory
__EXE_PATH=${__RAW_DATA_PATH_DIR} 
########################################################################
__INPUT_SAMPLE_DIR_List=(
#gene_promoterGenebody_UP_D_EX_10k_iv
memoryCD8_foxp3_common_peaks
Treg_foxp3_common_peaks
WT_CD4_foxp3_common_peaks
WT_CD8_foxp3_common_peaks
memoryCD8_solo_peaks
Treg_solo_peaks
WT_CD4_solo_peaks
WT_CD8_solo_peaks
#Overlap_Treg_vs_foxp3_13    
#Treg_solo_peaks_13
#Overlap_WT_CD8_vs_foxp3_13  
#WT_CD8_solo_peaks_13
)

__INPUT_SAMPLE_DIR_List=(
Sample_Treg_TCF1_20160827000
Sample_TKOCD4_TCF1_20160827000
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
RUN_CHOOSE
Alert_email=$?

echo ""
echo "__FASTQ_DIR_R1 __FASTQ_DIR_R2 are the READS_FULL_DIR FOR ANALYSIS"

for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	#RUN_Wig2BigWig ${__INPUT_SAMPLE_DIR_List[0]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[1]}
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[1]}
	#RUN_bed_intersect ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_CUT_Columns ${__INPUT_SAMPLE_DIR_List[i]} 1 3
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}_13
	
	
	break
done



	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}"
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

