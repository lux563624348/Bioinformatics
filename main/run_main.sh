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
__INPUT_PATH=/home/szhu/stat5/macs2
### Output DIRECTORY
__OUTPUT_PATH=/home/xli/Data_Processing/Haihui/ChIP_Seq_STAT5
__INPUT_SAMPLE_SET=(
dko_0h_1
dko_0h_2
dko_24h_1
dko_24h_2
dko_24h_3
dko_24h_4
dko_3h_1
dko_3h_2
dko_3h_3
dko_3h_4
wt_0h_1
wt_0h_2
wt_24h_1
wt_24h_2
wt_24h_3
wt_24h_4
wt_3h_1
wt_3h_2
wt_3h_3
wt_3h_4
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
Project_Name='CD8_STAT5'
##....................................................................##
### Download Raw Data
#FUNC_Download "ftp://ftp.admerahealth.com/19092-06" "19092-06"

##....................................................................##
### PARALLEL OPERTATION
echo "Parallel Operation have started"
for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_SET[*]} - 1); i++ ))  ### Loop Operation [Ref.1]
do
	RUN_BedGraph2BigWig ${__INPUT_PATH} ${__INPUT_SAMPLE_SET[i]}_treat_pileup ${Project_Name} ${SPECIES} ${Data_Provider} & pid=$!
	
	#RUN_SICER ${__INPUT_SAMPLE_SET[0]} ${__INPUT_SAMPLE_SET[1]} 200 ${Project_Name} ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Motif_Homer ${__INPUT_SAMPLE_SET[i]} 'null' ${SPECIES} ~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/RNA_seq/CuffDiff_Jun2018/Cuffdiff_Results/DKO_0h_vs_WT_0h/gene_exp.diff 'no' & pid=$!
	#REMOVE_REDUNDANCY_PICARD ${__INPUT_SAMPLE_SET[i]} & pid=$!
	#RUN_SRA2FASTQ ${__INPUT_SAMPLE_SET[i]}
	#PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fastq.gz' 'Pair'
	#RUN_HiC_Iterative_Mapping ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Trim_Galore_QC & pid=$!
	#RUN_FAST_QC & pid=$!
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
	#RUN_RPKM ${__INPUT_SAMPLE_SET[i]} bed 'customeized' 'yes' & pid=$!
	#RUN_MACS2 ${__INPUT_SAMPLE_SET[1]} 'Null' ${Project_Name} ${SPECIES} ${Data_Provider} 'bed' & pid=$!
	#RUN_TOPHAT ${__INPUT_SAMPLE_SET[i]} "TEST" ${SPECIES} ${Data_Provider} & pid=$!
	#RUN_Venn_Diagram ${__INPUT_PATH} 'bed'
	#RUN_ROSE_SUPER_Enhancer ${__INPUT_SAMPLE_SET[i]} ${SPECIES} 'customeized' & pid=$!
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_bed2fastq ${__INPUT_SAMPLE_SET[i]} ${SPECIES} & pid=$!
	#RUN_Reads_Profile "GeneBody" ${__INPUT_SAMPLE_Set[i]} ${SPECIES}

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
	PRE_READS_DIR ${__INPUT_SAMPLE_SET[i]} 'fq.gz' 'Pair'
	RUN_BOWTIE2 ${__INPUT_SAMPLE_SET[i]} ${SPECIES} ${Project_Name} ${Data_Provider} 'no' & pid=$!
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
