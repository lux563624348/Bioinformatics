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
__RAW_DATA_PATH_DIR=~/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/DNase_seq
#### Execution or Output directory
__EXE_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq
########################################################################

###Pool 1
__INPUT_SAMPLE_DIR_List=(
Sample_IgGForTLE3_20180702000
Sample_TCF1-TKO-20170610_20180702000
Sample_TCF1-Treg-20170610_20180702000
Sample_TCF1-nT-20170610_20180702000
Sample_TCF1-naiveCD8_20180702000
Sample_TCF1-sti72h_20180702000
Sample_TCFI-TKO_20180702000
Sample_TLE3-1_20180702000
Sample_TLE3-2_20180702000
Sample_TLE3-3_20180702000
Sample_TLE3-4_20180702000
fastq_screen_results
)
### Pool 2
__INPUT_SAMPLE_DIR_List=(
Sample_CD8-72hrs-sti_20180702000
Sample_CD8-TKO_20180702000
Sample_DKO-stat5b3h_20180702000
Sample_DKO-stat5b24h_20180702000
Sample_Ezh2-EKOCD4_20180702000
Sample_Ezh2-TLKO-Tfh_20180702000
Sample_Ezh2-ctrl-Tfh_20180702000
Sample_IgGforStat5b_20180702000
Sample_ctrl-stat5b3h-400_20180702000
Sample_ctrl-stat5b3h_20180702000 
Sample_ctrl-stat5b24h-400_20180702000
Sample_ctrl-stat5b24h_20180702000
fastq_screen_results
)


### Pool 3
__INPUT_SAMPLE_DIR_List=(
Sample_CD8-Naive_20180702000
Sample_K27Ac-TLE134-CD84_20180702000
Sample_K27ac-Gai_20180702000
Sample_K27ac-Zhao_20180702000
Sample_K27ac-shao_20180702000
Sample_K27me3-EKO-Tfh_20180702000
Sample_K27me3-ctrl-Tfh_20180702000
Sample_WT-input_20180702000
Sample_input-EKO-Tfh_20180702000
Sample_input-TLE134-CD84_20180702000 
Sample_input-ctrl-Tfh_20180702000
fastq_screen_results
)

__INPUT_SAMPLE_DIR_List=(
Intersection_Foxp3_5k_Treg_TCF1
Only_Foxp3_5k
Only_Treg_TCF1
)


__INPUT_SAMPLE_DIR_List=(
10_New_H3K9me3_WT    
11_New_H3K9me3_mir34bc_KO
5_H3K4me3_WT
7_H3K4me3_mir34bc_KO
8_H3K9me3_WT_April
8_H3K9me3_WT_Dec
9_H3K9me3_mir34bc_KO_April
9_H3K9me3_mir34bc_KO_Dec
)

# Pool 4
__INPUT_SAMPLE_DIR_List=(
Sample_WT-na1_20180709000
Sample_WT-na2_20180709000
Sample_WT-s1_20180709000
Sample_WT-s2_20180709000
Sample_dKO-na1_20180709000
Sample_dKO-na2_20180709000
Sample_dKO-s1_20180709000
Sample_dKO-s2_20180709000
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


### Parallel TEST
parallel_process(){
	break
	}

#for INPUT in ${__INPUT_SAMPLE_DIR_List[*]}; do RUN_Reads_Profile_Promoter 'TSS' ${INPUT} & done


for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Reads_Profile_Promoter 'TSS' ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]}
	FUNC_Download ${__INPUT_SAMPLE_DIR_List[i]}
	PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	RUN_FAST_QC
	#RUN_Venn_Diagram ${__EXE_PATH} 'bed'
	RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} "mm10"
	#RUN_RPKM ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_Reads_Profile_Promoter_genebody ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg" "mm10" "WT_Online_Ref"
	#RUN_BED2WIG ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	
	#RUN_Wig2BigWig ${__RAW_DATA_PATH_DIR} ${__INPUT_SAMPLE_DIR_List[i]} 'Tcf1' ${SPECIES} ${Data_Provider}
	
	#RUN_CUFFDIFF ${__INPUT_SAMPLE_DIR_List[*]}
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
	
	#### Parallel TEST
done
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[2]}
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[1]} ${__INPUT_SAMPLE_DIR_List[2]}
## SECOND LOOP
	#RUN_Reads_Profile "TSS"
	
	ps j > PPID.log

	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

