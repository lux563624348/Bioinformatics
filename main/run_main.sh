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
__RAW_DATA_PATH_DIR=~/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/DNase_seq/Raw_Reads
#### Execution or Output directory
__EXE_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNaseq_seq_RNA_Seq/Only_R1_TEST
########################################################################

#echo "INPUT_SAMPLE_DIR_List= (${__INPUT_SAMPLE_DIR_List[*]})"

########################################################################

########################################################################
##	MAIN BODY
########################################################################


__INPUT_SAMPLE_DIR_List=(
WT-na1_20180709000
WT-na2_20180709000
WT-s1_20180709000
WT-s2_20180709000
dKO-na1_20180709000
dKO-na2_20180709000
dKO-s1_20180709000
dKO-s2_20180709000
)

__INPUT_SAMPLE_BARCODE_List=(
ACAGTG
CAGATC
)

main() {
#### Saving DIR Check and Create
DIR_CHECK_CREATE $__RAW_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH

#### Email Alert
FUNC_CHOOSE_EMAIL_ALERT
Alert_email=$?

echo ""
echo "__FASTQ_DIR_R1 __FASTQ_DIR_R2 are the READS_FULL_DIR FOR ANALYSIS"

###
SPECIES='mm9'
Data_Provider='Haihui'
####

for (( i = 0; i <= $(expr ${#__INPUT_SAMPLE_DIR_List[*]} - 1); i++ ))
do
	PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} 'fastq.gz' 'R1'
	#RUN_FAST_QC
	#RUN_HomerTools 'Restriction_Enzyme' ${__INPUT_SAMPLE_BARCODE_List[i]} &
	RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES} "CD8-HP_R1_Test" ${Data_Provider} &
	
	#break
	#RUN_HomerTools 'Restriction_Enzyme' ${__INPUT_SAMPLE_BARCODE_List[i]}
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES} "CD8-HP-HiC_seq_bowtie2_Only_R2" ${Data_Provider} 
	#break
	#RUN_RPKM ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	#RUN_CELLRANGER ${__INPUT_SAMPLE_DIR_List[i]} "Hdac" "mm10"
	#RUN_Bed2BigBed ${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[i]} ${__INPUT_SAMPLE_DIR_List[i]} "CD8-HP-Only_R1" ${SPECIES} ${Data_Provider}
	#RUN_BedGraph2BigWig ${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[i]} ${__INPUT_SAMPLE_DIR_List[i]} "Bed2bdg2bigwig" ${SPECIES} ${Data_Provider}
	#RUN_Island_Filtered_Reads ${__INPUT_SAMPLE_DIR_List[i]} 'bedpe' &
	#RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[i]} 'Null' 'TCF1_ChIP_Seq_Macs2_SPMR' ${SPECIES} ${Data_Provider} 'bam' &
	#FUNC_Download ${__INPUT_SAMPLE_DIR_List[i]}

	#RUN_FAST_QC &
	
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Reads_Profile_Promoter 'TSS' ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg_RNA-seq_201806" "mm9" "Haihui"
	#RUN_Venn_Diagram ${__EXE_PATH} 'bed'
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} "mm10"
	
	#RUN_Reads_Profile_Promoter_genebody ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg" "mm10" "WT_Online_Ref"
	#RUN_BED2WIG ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	
	#RUN_Wig2BigWig ${__RAW_DATA_PATH_DIR}/${__INPUT_SAMPLE_DIR_List[i]} ${__INPUT_SAMPLE_DIR_List[i]} 'CD8-HP-DNase_seq' ${SPECIES} ${Data_Provider}
	
	#RUN_CUFFDIFF ${__INPUT_SAMPLE_DIR_List[*]}
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

