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
__RAW_DATA_PATH_DIR=/home/lxiang/Raw_Data/Paul/34bc/raw_data/8_H3K9me3_WT_April/
__RAW_Files_Name=(H3K9me3_WT_R1_goodreads.fastq H3K9me3_WT_R2_goodreads.fastq)

__Target_DATA_PATH_DIR=/home/lxiang/Data/Paul/34bc/MMLVless_data_protocol_B/8_H3K9me3_WT_April/bowtie2_results/
__Target_Files_Name=(conc_aligned_R1.fastq conc_aligned_R2.fastq)

__EXE_PATH=/home/lxiang/Data/Paul/34bc/Raw_data_Directly_remove_MMLV/


#### TEST
#__RAW_DATA_PATH_DIR=/home/lxiang/Data/Paul/34bc/Raw_data_Directly_remove_MMLV/TEST/
#__RAW_Files_Name=(pMXs-Sox2-full_R1.fastq)
#__Target_DATA_PATH_DIR=${__RAW_DATA_PATH_DIR}
#__Target_Files_Name=(target.fastq)



########################################################################
#echo "__INPUT_SAMPLE_DIR_List=(1_input_Bruce4 2_H3K4me3_Bruce4 3_H3K9me3_Bruce4 4_input_WT 5_H3K4me3_WT 6_input_mir34bc_KO 7_H3K4me3_mir34bc_KO 8_H3K9me3_WT_Dec 8_H3K9me3_WT_April 9_H3K9me3_mir34bc_KO_Dec 9_H3K9me3_mir34bc_KO_April 10_New_H3K9me3_WT 11_New_H3K9me3_mir34bc_KO)"

SAMPLE_NUM=${#__RAW_Files_Name[*]}
########################################################################

########################################################################
##	MAIN BODY
########################################################################

main() {
#### Saving DIR Check and Create
DIR_CHECK_CREATE $__RAW_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH

echo ""
echo ""
########################################################################

#### Raw_data
rows_raw_data=$(cat ${__RAW_DATA_PATH_DIR}${__RAW_Files_Name[0]} | wc -l)
reads_raw_data=$(expr $rows_raw_data / 4)
echo "Number of raw reads: $reads_raw_data"
input_data_path=${__RAW_DATA_PATH_DIR}${__RAW_Files_Name[0]}
echo 'wc -l ${input_data_path}'


#### Target_data
rows_target_data=$(cat ${__Target_DATA_PATH_DIR}${__Target_Files_Name[0]} | wc -l)
reads_target_data=$(expr $rows_target_data / 4)
echo "Number of target reads : $reads_target_data"
input_target_path=${__Target_DATA_PATH_DIR}${__Target_Files_Name[0]}
echo 'wc -l ${input_target_path}'

######
i=0
while [ $i -lt ${reads_raw_data} ]
#for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	rows_start=$(expr 4 \* $i + 1)
	rows_end=$(expr $rows_start + 3)
	Query_name=$(sed -n "${rows_start}p" ${input_data_path})
	
	j=0
	while  [ $j -lt ${reads_target_data} ]
	do
	target_rows_start=$(expr 4 \* $j + 1)
	target_rows_end=$(expr ${target_rows_start} + 3)
	target_Query_name=$(sed -n "${target_rows_start}p" ${input_target_path})
	
	#### Comparing part
	if [ ${Query_name:20 :20} = ${target_Query_name:20 :20} ];then
	echo "Skip Query name: ${Query_name}"
	break
	fi
	
	if [ $(expr $j + 1 ) -eq ${reads_target_data} ];then
	sed -n "${rows_start},${rows_end}p" ${input_data_path} >> Output_R1.fastq
	fi
	
	j=$(expr $j + 1)
	done
	
	i=$(expr $i + 1)
done
echo "End Date: `date`"
echo -e "\a FINISHED ALERT !"
#EMAIL_ME $Start_Date
}

##### Following Line is very IMPORTANT 
main "$@"

