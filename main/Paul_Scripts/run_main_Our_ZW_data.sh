#!/bin/bash
#set -e	#### Terminating the script if any command exited with a nonzero exit status
set -u	#### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.
### functions.sh contains all modules and functions necessary in this pipeline, has to be together with run_main.sh
source ./functions.sh
export PATH="/public/slst/home/fanyh/anaconda3/envs/py3_lx/bin:$PATH"
## USAGE:
#	nohup bash run_main.sh "Process_NAME" > out.log &
########################################################################
## 02/22/2019
## By Xiang Li,
## lux@gwu.edu
## Li's Lab
## Version.beta
########################################################################

########################################################################
  #  means NEED MODIFICATION. "VERY IMPORTANT INFORMATION"
  ## means Title_level 1.
 ### means Title_level 2.
#### means comment.
########################################################################

## Pickup a NAME for this one time pipeline.
Process_NAME="XL_pipeline"
echo "-----------------------------------------------------------------"

echo "Import functions.sh Completed!"
echo "Your Pipeline Named ${Process_NAME} start at `date`"
Start_Date=`date`
echo "-----------------------------------------------------------------"
#			NEED MODIFICATION FOR DIFFERENT PROJECT
########################################################################
## GLOBAL VARIABLES
########################################################################
### INPUT DIRECTORY
__INPUT_PATH=/slst/home/fanyh/Raw_Data/Totipotent/RNAseq
### Output DIRECTORY
__OUTPUT_PATH=/slst/home/fanyh/Data/JGA
__INPUT_SAMPLE_SET=(
EVT-1
hTSC-1
hESC-cas9-TB-1
EVT-D3-NC-1
EVT-D3-F3-2
EVT-D3-F2-1
EVT-D3-F1-2
hESC-cas9-2
hESC-cas9-1
hTSC-2
hESC-cas9-TB-2
EVT-D3-NC-2_FRA
EVT-D3-F3-1_FRA
EVT-D3-F2-2_FRA
EVT-D3-F1-1_FRA
EVT-2_FRAS20241
)
#### Saving DIR Check and Create
DIR_CHECK_CREATE ${__OUTPUT_PATH} ${__INPUT_PATH}
########################################################################
#### Function of Email Alert, yes or no?
#FUNC_CHOOSE_EMAIL_ALERT
#Alert_email=$?
##	MAIN BODY
########################################################################
main() {
echo "-----------------------------------------------------------------"
echo "$(date "+%Y-%m-%d %H:%M") Start Processing....."
##....................................................................##
### Key Parameters
SPECIES='hg38'
Data_Provider='Paul'
Project_Name='Totipotent'
##....................................................................##
### Download Raw Data
#FUNC_Download "ftp://ftp.admerahealth.com/19092-06" "19092-06"

##....................................................................##
### PARALLEL OPERTATION
echo "Parallel Operation have started"
#conda activate py3_lx
for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	PRE_READS_DIR ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'SRA'
	RUN_FAST_QC & pid=$!
	RUN_HISAT2 ${__INPUT_SAMPLE_SET[i]} ${Project_Name} ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Trim_Galore_QC & pid=$!
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	#RUN_MACS2 ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'Null' ${Project_Name} ${SPECIES} ${Data_Provider} 'bampe' & pid=$!
	#RUN_ROSE_SUPER_Enhancer ${__INPUT_SAMPLE_SET[i]} ${SPECIES} 'customeized' & pid=$!
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_bed2fastq ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Reads_Profile "GeneBody" ${__INPUT_SAMPLE_Set[i]} ${SPECIES}
	#REMOVE_REDUNDANCY_PICARD ${__INPUT_SAMPLE_SET[i]} & pid=$!
	
	
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
	break
	PRE_READS_DIR ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'Pair'
	#RUN_Quant_IRI ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} & pid=$!
	#PRE_READS_DIR ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'Pair'
	RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	#RUN_MACS2 ${__INPUT_SAMPLE_SET[i]} 'Null' ${Project_Name} ${SPECIES} ${Data_Provider} 'bampe' & pid=$!
	PID_LIST+=" $pid";
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
	
#	if [ ${Alert_email} == 0 ];then
#	EMAIL_ME "${Start_Date}" ${Process_NAME}
#	fi
}
##### Following Line is very IMPORTANT 
main "$@"

#[1] https://stackoverflow.com/questions/10909685/run-parallel-multiple-commands-at-once-in-the-same-terminal
#[2]
#[3]  
