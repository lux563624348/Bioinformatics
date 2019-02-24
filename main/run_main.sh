#!/bin/bash
set -e	#### Terminating the script if any command exited with a nonzero exit status
set -u	#### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.

## USAGE:
#	nohup bash run_main.sh "Process_NAME" > out.log &
########################################################################
## 02/22/2019
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

## Pickup a NAME for this one time pipeline.
Process_NAME=${1}
echo "-----------------------------------------------------------------"
### functions.sh contains all modules and functions necessary in this pipeline, has to be together with run_main.sh
source ./functions.sh
echo "Import functions.sh Completed!"
echo "Your Pipeline Named ${Process_NAME} start at `date`"
Start_Date=`date`
echo "-----------------------------------------------------------------"
#			NEED MODIFICATION FOR DIFFERENT PROJECT
########################################################################
## GLOBAL VARIABLES
########################################################################
### Output DIRECTORY
__OUTPUT_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq/MACS2_Results/bampe/Union_all/version.beta
### INPUT DIRECTORY
__INPUT_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq/MACS2_Results/bampe/Union_all
#### Unique Keyword of the library you want to put into pipeline.
__INPUT_SAMPLE_SET=(
12532_TCF1+_TCF1_MOTIF_-_Union_peak_For_HTseq
1297_Tcf1_Motif_+_Union_peak_For_HTseq
)
#### Saving DIR Check and Create
DIR_CHECK_CREATE ${__OUTPUT_PATH} ${__INPUT_PATH}
########################################################################
#### Function of Email Alert, yes or no?
FUNC_CHOOSE_EMAIL_ALERT
Alert_email=$?
##	MAIN BODY
########################################################################
main() {
echo "-----------------------------------------------------------------"
echo "$(date "+%Y-%m-%d %H:%M") Start Processing....."
##....................................................................##
### Key Parameters
SPECIES='mm9'
Data_Provider='Haihui'
##....................................................................##
### Download Raw Data
#FUNC_Download "http://shang.phys.gwu.edu/Test_Data" "test_Pipeline"

##....................................................................##
### PARALLEL OPERTATION
echo "Parallel Operation have started";
for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	#PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fastq.gz' 'Pairs'
	#RUN_FAST_QC &
	#RUN_TOPHAT ${__INPUT_SAMPLE_SET[i]} "TEST" ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Venn_Diagram ${__INPUT_PATH} 'bed'
	#RUN_ROSE_SUPER_Enhancer ${__INPUT_SAMPLE_Set[0]}
	RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Reads_Profile "GeneBody" ${__INPUT_SAMPLE_Set[i]} ${SPECIES}
	#RUN_SICER ${__INPUT_SAMPLE_List[1]} ${__INPUT_SAMPLE_Set[0]} 400 "CD8-K27Ac" ${SPECIES} ${Data_Provider} & pid=$!
	PID_LIST+=" $pid";
	#break
	#break
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
done
echo "wait ${PID_LIST}}....................................."
wait ${PID_LIST}
echo "Parallel Operation have finished";
##....................................................................##
### Single Operation
	#RUN_CELLRANGER ${__INPUT_SAMPLE_DIR_List[15]} "Hdac" "mm10"
echo "Single Operation have finished";

##....................................................................##
	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}
##### Following Line is very IMPORTANT 
main "$@"

#[1] https://stackoverflow.com/questions/10909685/run-parallel-multiple-commands-at-once-in-the-same-terminal
#[2]
#[3]  
