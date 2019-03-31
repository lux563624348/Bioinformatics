#!/bin/bash
#set -e	#### Terminating the script if any command exited with a nonzero exit status
set -u	#### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.
### functions.sh contains all modules and functions necessary in this pipeline, has to be together with run_main.sh
source ./functions.sh

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

echo "Import functions.sh Completed!"
echo "Your Pipeline Named ${Process_NAME} start at `date`"
Start_Date=`date`
echo "-----------------------------------------------------------------"
#			NEED MODIFICATION FOR DIFFERENT PROJECT
########################################################################
## GLOBAL VARIABLES
########################################################################
### Output DIRECTORY
__OUTPUT_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/1903_Haihui_link1
### INPUT DIRECTORY
__INPUT_PATH=~/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/1903_Haihui_link1
#### Unique Keyword of the library you want to put into pipeline.
__INPUT_SAMPLE_SET=(
Input_20190301000
naive-WT-CD8-1_20190301000
naive-WT-CD8-2_20190301000 
naive-dKO-CD8-2_20190301000
stim-WT-CD8-1_20190301000
stim-WT-CD8-2_20190301000
stim-dKO-CD8-2_20190301000
ZC-less-stringent-wash_20190301000
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
Project_Name='CD8-HP'
##....................................................................##
### Download Raw Data
#FUNC_Download "http://dnacore454.healthcare.uiowa.edu/20190318-0351_Xue_QiangQS-C-PoolrkaxYXeCWDmgtBFhrIbqeRTPeHsZPJpbbHIKkDRO/results/" "1903_Haihui_link3"

##....................................................................##
### PARALLEL OPERTATION
echo "Parallel Operation have started"
for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	#RUN_Motif_Homer ${__INPUT_SAMPLE_SET[0]} ${__INPUT_SAMPLE_SET[2]} ${SPECIES} "No Expression" 'no' & pid=$!
	#RUN_Motif_Homer ${__INPUT_SAMPLE_SET[1]} ${__INPUT_SAMPLE_SET[3]} ${SPECIES} "No Expression" 'no' & pid=$!
	#RUN_SRA2FASTQ ${__INPUT_SAMPLE_SET[i]}
	PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fastq.gz' 'Pair'
	break
	RUN_Trim_Galore_QC & pid=$!
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	#RUN_MACS2 ${__INPUT_SAMPLE_SET[i]} "Null" WT_CFCT ${SPECIES} ${Data_Provider} 'bam' & pid=$!
	#RUN_TOPHAT ${__INPUT_SAMPLE_SET[i]} "TEST" ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Venn_Diagram ${__INPUT_PATH} 'bed'
	#RUN_ROSE_SUPER_Enhancer ${__INPUT_SAMPLE_Set[0]}
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Reads_Profile "GeneBody" ${__INPUT_SAMPLE_Set[i]} ${SPECIES}
	#RUN_SICER ${__INPUT_SAMPLE_List[1]} ${__INPUT_SAMPLE_Set[0]} 400 "CD8-K27Ac" ${SPECIES} ${Data_Provider} & pid=$!
	PID_LIST+=" $pid";
	#break
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
done
echo "wait ${PID_LIST}....................................."
wait ${PID_LIST}
echo "Parallel Operation have finished";

for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'Pair'
	RUN_MACS2 ${__INPUT_SAMPLE_SET[i]} ${__INPUT_SAMPLE_SET[i]} "xx" ${SPECIES} ${Data_Provider} 'bam' & pid=$!
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	PID_LIST+=" $pid";
	#break
done
echo "wait ${PID_LIST}....................................."
wait ${PID_LIST}
echo "Parallel Operation have finished";

##....................................................................##
### Single Operation
	#RUN_CELLRANGER ${__INPUT_SAMPLE_DIR_List[15]} "Hdac" "mm10"
#echo "Single Operation have finished";

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
