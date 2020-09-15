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
### INPUT DIRECTORY
__INPUT_PATH=/home/xli/Data/Haihui/mm10_Xin
### Output DIRECTORY
__OUTPUT_PATH=/home/xli/Data/Haihui/mm10_Xin
__INPUT_SAMPLE_SET=(
#naiveCD8_TCF1 ##
SRR5217327 ##GSM2471915   >  Naïve CD8 CBF ChIP-Seq
SRR5217328 #GSM2471916    >  CBFβ-deficient naive CD8 CBF ChIP-Seq
SRR955620 #GSM1214544     >  Runx3_CD8+T_CHIP_Runx3_IP1 file 1
SRR955621 #GSM1214545     >  Runx3_CD8+T_CHIP_Runx3_IP1 file 2
SRR955622 #GSM1214546     >  Runx3_CD8+T_CHIP_Runx3_IP2
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
SPECIES='mm10'
Data_Provider='Haihui'
Project_Name='Xin_mm10'
##....................................................................##
### Download Raw Data
#FUNC_Download "ftp://ftp.admerahealth.com/19092-06" "19092-06"

##....................................................................##
### PARALLEL OPERTATION
echo "Parallel Operation have started"
#conda activate py3_lx
for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	#RUN_BedGraph2BigWig ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]}_treat_pileup ${Project_Name} ${SPECIES} ${Data_Provider} & pid=$!
	
	#RUN_SICER ${__INPUT_SAMPLE_SET[0]} ${__INPUT_SAMPLE_SET[1]} 200 ${Project_Name} ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Motif_Homer ${__INPUT_SAMPLE_SET[i]} 'Null' ${SPECIES} 'No Expression' 'yes' & pid=$!
	#RUN_Island_Filtered_Reads ${__INPUT_SAMPLE_SET[i]} 'bed' & pid=$!
	#break
	#REMOVE_REDUNDANCY_PICARD ${__INPUT_SAMPLE_SET[i]} & pid=$!
	#RUN_SRA2FASTQ ${__INPUT_SAMPLE_SET[i]} & pid=$!
	PRE_READS_DIR ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'Pair'
	#PRE_READS_DIR ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'fastq.gz' 'Pair'
	#RUN_FAST_QC & pid=$!
	#RUN_HiC_Iterative_Mapping ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'SRA'
	#RUN_HiC_Iterative_Mapping ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Trim_Galore_QC & pid=$!
	#RUN_FAST_QC & pid=$!
	RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	#RUN_RPKM ${__INPUT_SAMPLE_SET[i]} bam 'customeized' 'yes' & pid=$!
	#RUN_MACS2 ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]} 'Null' ${Project_Name} ${SPECIES} ${Data_Provider} 'bampe' & pid=$!
	#RUN_TOPHAT ${__INPUT_SAMPLE_SET[i]} "szhu" ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Venn_Diagram ${__INPUT_PATH} 'bed'
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
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}
##### Following Line is very IMPORTANT 
main "$@"

#[1] https://stackoverflow.com/questions/10909685/run-parallel-multiple-commands-at-once-in-the-same-terminal
#[2]
#[3]  
