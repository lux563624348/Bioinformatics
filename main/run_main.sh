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
echo "-----------------------------------------------------------------"
echo "Start Date: `date`"
echo "-----------------------------------------------------------------"
Start_Date=`date`

Process_NAME=${1}
source ./functions.sh
echo "Import functions.sh Completed!" 
#			NEED MODIFICATION FOR DIFFERENT PROJECT

########################################################################
## GLOBAL VARIABLES
########################################################################
__RAW_DATA_PATH_DIR=~/cloud_research/PengGroup/XLi/Raw_Data/TEST_Raw_Data/RNA_seq
#### Execution or Output directory
__EXE_PATH=~/cloud_research/PengGroup/XLi/Raw_Data/TEST_Data/RNA_seq
#__EXE_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq
########################################################################
########################################################################
##	MAIN BODY
########################################################################
__INPUT_SAMPLE_List=(
Sample_FY1_20190114000
Sample_FY2_20190114000
Sample_FY3_20190114000
Sample_FY4_20190114000
Sample_FY5_20190114000
Sample_SP1_20190114000
Sample_SP2_20190114000
Sample_SP3_20190114000
Sample_SP4_20190114000
Sample_SP5_20190114000
Sample_SP6_20190114000
Sample_SP7_20190114000
Sample_SP8_20190114000
Sample_SP9_20190114000
Sample_SP10_20190114000
)

__INPUT_SAMPLE_List=(
Sample_GK1_20190114000
Sample_GK2_20190114000
Sample_GK3_20190114000
Sample_GK4_20190114000
Sample_GK5_20190114000
Sample_GK6_20190114000
Sample_GK7_20190114000
Sample_GK8_20190114000
Sample_GK9_20190114000
Sample_GK10_20190114000
Sample_GK11_20190114000
)

__INPUT_SAMPLE_List=(
FY1_test
)

main() {
#### Saving DIR Check and Create
DIR_CHECK_CREATE $__RAW_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH

#### Email Alert
FUNC_CHOOSE_EMAIL_ALERT
Alert_email=$?

echo "-----------------------------------------------------------------"
echo "$(date "+%Y-%m-%d %H:%M") Start Processing....."

###
SPECIES='mm10'
Data_Provider='Haihui'
####

for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_List[*]} - 1); i++ ))
do
	FUNC_Download ${__INPUT_SAMPLE_List[i]} "http://shang.phys.gwu.edu/Test_Data"
	#RUN_SRA2FASTQ ${__INPUT_SAMPLE_DIR_List[i]} 
	PRE_READS_DIR ${__INPUT_SAMPLE_List[i]} 'fastq.gz' 'Pairs'
	RUN_FAST_QC  
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_List[i]} ${SPECIES} "Pre_Tfh_Th1" ${Data_Provider} 'no' &
	RUN_TOPHAT ${__INPUT_SAMPLE_List[i]} "TEST" ${SPECIES} ${Data_Provider} 
	#RUN_CUFFDIFF ${__INPUT_SAMPLE_DIR_List[*]}
	#PRE_READS_DIR ${__INPUT_SAMPLE_List[0]} 'fastq.gz' 'Pairs'
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_List[0]} ${SPECIES} "CD8_hm_ChIPseq" ${Data_Provider} 'no'  
	#RUN_MACS2 ${__INPUT_SAMPLE_List[0]} ${__INPUT_SAMPLE_DIR_List[1]} 'CD8_hm_ChIPseq' ${SPECIES} ${Data_Provider} 'bam'
	
	#RUN_ROSE_SUPER_Enhancer ${__INPUT_SAMPLE_List[i]} ${SPECIES}
	
	#RUN_Motif_Homer ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[1]} ${SPECIES} 'yes'

	#RUN_RPKM ${__INPUT_SAMPLE_DIR_List[1]} 'bed' ${SPECIES} 
	#break
	
	#RUN_MACS2_Diff ${__INPUT_SAMPLE_DIR_List[3]} ${__INPUT_SAMPLE_DIR_List[1]}
	#RUN_CELLRANGER ${__INPUT_SAMPLE_DIR_List[i]} "Hdac" "mm10"
	#RUN_Island_Filtered_Reads ${__INPUT_SAMPLE_DIR_List[i]} 'bedpe' &
	#RUN_Reads_Profile "GeneBody" ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES} &
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[i]} 'Null' 'CD8-HP_DNaseq_MACS2_Merge_Replicates' ${SPECIES} ${Data_Provider} 'bampe' &
	#RUN_bed2fastq ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}

	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Reads_Profile_Promoter 'TSS' ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg_RNA-seq_201806" "mm9" "Haihui"
	#RUN_Venn_Diagram ${__RAW_DATA_PATH_DIR} 'bed'
	#break
	
	#RUN_Reads_Profile_Promoter_genebody ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg" "mm10" "WT_Online_Ref"
	#RUN_BED2WIG ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	
	#RUN_Wig2BigWig ${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[i]} ${__INPUT_SAMPLE_DIR_List[i]} 'CD8-HP-DNase_seq' ${SPECIES} ${Data_Provider}
	
	#
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
	
done
	
	#RUN_CELLRANGER ${__INPUT_SAMPLE_DIR_List[15]} "Hdac" "mm10"
	
	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

