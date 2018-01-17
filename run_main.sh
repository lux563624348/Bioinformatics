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
__RAW_DATA_PATH_DIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/Project_Xue_TCF1_620/MACS2_Diff_results/naiveCD8_TCF1_vs_MemoryCD8_TCF1_FDR_E-5
#### Execution or Output directory
__EXE_PATH=${__RAW_DATA_PATH_DIR}
#__EXE_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/DN3
########################################################################
echo "__INPUT_SAMPLE_DIR_List=(1_input_Bruce4 2_H3K4me3_Bruce4 3_H3K9me3_Bruce4 4_input_WT 5_H3K4me3_WT 6_input_mir34bc_KO 7_H3K4me3_mir34bc_KO 8_H3K9me3_WT_Dec 8_H3K9me3_WT_April 9_H3K9me3_mir34bc_KO_Dec 9_H3K9me3_mir34bc_KO_April 10_New_H3K9me3_WT 11_New_H3K9me3_mir34bc_KO)"

__INPUT_SAMPLE_DIR_List=(
Sample_DN3_LEF1_20160901000  #DN3 thymocytes
Sample_Thy_LEF1_20160901000  #whole thymocytes
Sample_vavTKO_LEF1_20160901000 # TKO whole thymocytes
Sample_LKO_LFF1_20160901000 # LKO whole thymocytes
)

__INPUT_SAMPLE_DIR_List=(
Sample_naiveCD8_TCF1_20160827000
Sample_memoryCD8_TCF1_20160827000
Sample_TKOCD8_TCF1_20160827000
)

__INPUT_SAMPLE_DIR_List=(
Sample_naiveCD8_TCF1_20160827000_island_reads
Sample_memoryCD8_TCF1_20160827000_island_reads
naiveCD8_CONTRO_island_reads
MemoryCD8_CONTRO_island_reads
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
	RUN_MACS2_Diff "naiveCD8_TCF1" ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[2]} \
	"MemoryCD8_TCF1" ${__INPUT_SAMPLE_DIR_List[1]} ${__INPUT_SAMPLE_DIR_List[3]} 0.00001
	break
	#PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]}
	
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset __FASTQ_DIR_R1 __FASTQ_DIR_R2
done


## SECOND LOOP
for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[1]}
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

