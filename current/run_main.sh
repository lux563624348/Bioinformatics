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
#__RAW_DATA_PATH_DIR=$(pwd)
__RAW_DATA_PATH_DIR=/home/lxiang/cloud_research/PengGroup/XLi/Raw_Data/Haihui/Tcf1/Project_Xue_TCF1_620
#### Execution or Output directory
__EXE_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620
########################################################################


__INPUT_SAMPLE_DIR_List=(
Sample_CD4_TCF1_20160827000 #0
Sample_Treg_TCF1_20160827000   #1
Sample_TKOCD4_TCF1_20160827000   #2
Sample_naiveCD8_TCF1_20160827000   #3
Sample_memoryCD8_TCF1_20160827000
Sample_TKOCD8_TCF1_20160827000    #5
Sample_DN3_TCF1_20160827000 #6
Sample_vavTKO_TCF1_20160827000 #7
Sample_CD4-IgG_20180212000  #8
Sample_CD4-TCF1_20180212000 #9
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
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]}
	
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[2]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[1]} ${__INPUT_SAMPLE_DIR_List[2]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[8]} ${__INPUT_SAMPLE_DIR_List[2]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[9]} ${__INPUT_SAMPLE_DIR_List[2]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[3]} ${__INPUT_SAMPLE_DIR_List[5]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[4]} ${__INPUT_SAMPLE_DIR_List[5]}
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[6]} ${__INPUT_SAMPLE_DIR_List[7]}
	
	
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

