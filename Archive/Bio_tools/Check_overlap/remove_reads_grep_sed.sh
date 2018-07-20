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
__RAW_DATA_PATH_DIR=/home/lxiang/Raw_Data/Paul/34bc/raw_data/9_H3K9me3_mir34bc_KO_Dec/
__RAW_Files_Name=(H3K9me3_mir34bc_KO_R1_goodreads.fastq H3K9me3_mir34bc_KO_R2_goodreads.fastq)

__Target_DATA_PATH_DIR=/home/lxiang/Data/Paul/34bc/MMLVless_data_protocol_B/9_H3K9me3_mir34bc_KO_Dec/bowtie2_results/
__Target_Files_Name=(conc_aligned_R1.fastq conc_aligned_R2.fastq)

__EXE_PATH=/home/lxiang/Data/Paul/34bc/Raw_data_Directly_remove_MMLV/9_H3K9me3_mir34bc_KO_Dec/


#### TEST
#__RAW_DATA_PATH_DIR=/home/lxiang/Data/Paul/34bc/Raw_data_Directly_remove_MMLV/TEST/
#__RAW_Files_Name=(pMXs-Sox2-full_R1.fastq)
#__Target_DATA_PATH_DIR=/home/lxiang/Data/Paul/34bc/Raw_data_Directly_remove_MMLV/TEST/
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
DIR_CHECK_CREATE $__Target_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH

echo ""
########################################################################


for ((i = 0; i <= $(expr $SAMPLE_NUM - 1); i++))
do
#### Raw_data
	rows_raw_data=$(cat ${__RAW_DATA_PATH_DIR}${__RAW_Files_Name[i]} | wc -l)
	reads_raw_data=$(expr $rows_raw_data / 4)
	echo "Number of raw reads: $reads_raw_data"
	input_data_path=${__RAW_DATA_PATH_DIR}${__RAW_Files_Name[i]}


#### Target_data
	rows_target_data=$(cat ${__Target_DATA_PATH_DIR}${__Target_Files_Name[i]} | wc -l)
	reads_target_data=$(expr $rows_target_data / 4)
	echo "Number of target reads : $reads_target_data"
	input_target_path=${__Target_DATA_PATH_DIR}${__Target_Files_Name[i]}

######
	Q_name_ref=@HS2
	
	Output_number=$(expr $i + 1 )
	Location_Out_name="Location.tem"
	Output_name="Output_R${Output_number}.fastq"
	
	grep --no-group-separator ${Q_name_ref} ${input_target_path} > Q_name_target.fastq
	grep --no-group-separator -xnf Q_name_target.fastq ${input_data_path} |cut -f1 -d:> ${Location_Out_name}
	
	rm Q_name_target.fastq
	
	rows_start=1
	for number in ${Location_Out_name}
	do
	if [ $number == 1 ];then
	rows_start=5
	continue
	else
	rows_end=$(expr $number - 1)
	if [ $rows_end -lt $rows_start ];then
	rows_start=$(expr ${number} + 4 )
	continue
	fi
	
	echo "Output: start: ${rows_start}, Ends: ${rows_end} "
	sed -n "${rows_start},${rows_end}p" ${input_data_path} >> ${Output_name}
	rows_start=$(expr ${number} + 4 )
	fi
	done
	
	Total_lines=$(cat $input_data_path | wc -l)
	echo "Output: start: ${rows_start}, Ends: ${Total_lines} "
	sed -n "${rows_start},${Total_lines}p" ${input_data_path} >> ${Output_name}
done


echo "End Date: `date`"
echo -e "\a FINISHED ALERT !"
#EMAIL_ME $Start_Date
}

##### Following Line is very IMPORTANT 
main "$@"

