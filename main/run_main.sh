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
Process_NAME=${1}


source ./functions.sh
echo "Import functions.sh" 

#			NEED MODIFICATION FOR DIFFERENT PROJECT

########################################################################
## GLOBAL VARIABLES
########################################################################
__RAW_DATA_PATH_DIR=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Hdac/Treg_HistoneMarks_ChIPSeq_Sep2016/bowtie_analysis
#### Execution or Output directory
__EXE_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Hdac/Treg_HistoneMarks_ChIPSeq_Sep2016/bowtie
########################################################################

__INPUT_SAMPLE_DIR_List=(
Hdac12_KO_Treg_D20_H3K27Ac
Hdac12_KO_Treg_D20_H3K9Ac
WT_Treg_D20_H3K27Ac
WT_Treg_D20_H3K9Ac
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
#echo "__FASTQ_DIR_R1 __FASTQ_DIR_R2 are the READS_FULL_DIR FOR ANALYSIS"

for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]}
	#FUNC_Download ${__INPUT_SAMPLE_DIR_List[i]}
	#PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	#RUN_FAST_QC
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} "mm9"
	RUN_RPKM ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Vlad" "mm9" "Haihui_Vlad"
	#RUN_BED2WIG ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	
	#SPECIES="mm9"
	#Data_Provider="Ref_GSM46662"
	#RUN_Wig2BigWig ${__RAW_DATA_PATH_DIR} ${__INPUT_SAMPLE_DIR_List[i]} 'Tcf1' ${SPECIES} ${Data_Provider}
	
	#RUN_CUFFDIFF ${__INPUT_SAMPLE_DIR_List[*]}
	
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
done

## SECOND LOOP
	#RUN_Reads_Profile "TSS"
#RUN_SICER

	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

