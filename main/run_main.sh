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
__RAW_DATA_PATH_DIR=~/cloud_research/PengGroup/XLi/Data/Haihui/Ezh2/ChIP_seq/Jul2018/Bowtie2_Results_filtered_S2/Bowtie2_Results
#### Execution or Output directory
__EXE_PATH=~/cloud_research/PengGroup/XLi/Data/Haihui/Ezh2/ChIP_seq/Jul2018
########################################################################

###Pool 1
__INPUT_SAMPLE_DIR_List=(
Sample_Ctrl-D1_20180801000
Sample_Ctrl-D1_20180801001
Sample_Ctrl-D1_20180801002
Sample_Ctrl-D1_20180801003
Sample_Ctrl-DO_20180801000
Sample_Ctrl-DO_20180801001
Sample_Ctrl-DO_20180801002
Sample_Ctrl-DO_20180801003
Sample_Dko-D1_20180801000
Sample_Dko-D1_20180801001
Sample_Dko-D1_20180801002
Sample_Dko-D1_20180801003
Sample_Dko-DO_20180801000
Sample_Dko-DO_20180801001
Sample_Dko-DO_20180801002
Sample_Dko-DO_20180801003
)

__INPUT_SAMPLE_DIR_List=(
EKO1_Tfh
EKO2_Tfh
WT1_Tfh
WT2_Tfh
)

__INPUT_SAMPLE_DIR_List=(
Sample_Ezh2-ctrl-Tfh_20180702000  #0
Sample_Ezh2-TLKO-Tfh_20180702000  #1
Sample_Ezh2-EKOCD4_20180702000    #2
Sample_input-ctrl-Tfh_20180702000 #3
Sample_input-EKO-Tfh_20180702000  #4
Sample_K27me3-ctrl-Tfh_20180702000 #5
Sample_K27me3-EKO-Tfh_20180702000 #6
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
echo "__FASTQ_DIR_R1 __FASTQ_DIR_R2 are the READS_FULL_DIR FOR ANALYSIS"


for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[2]} 'Ezh2_ChIP_seq_201807' &
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[1]} ${__INPUT_SAMPLE_DIR_List[2]} 'Ezh2_ChIP_seq_201807' &
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[3]} ${__INPUT_SAMPLE_DIR_List[5]} 'Ezh2_H3K27me3_201807' &
	RUN_MACS2 ${__INPUT_SAMPLE_DIR_List[4]} ${__INPUT_SAMPLE_DIR_List[6]} 'Ezh2_H3K27me3_201807'
	#FUNC_Download ${__INPUT_SAMPLE_DIR_List[i]}
	#PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} 'fastq.gz'
	#RUN_FAST_QC &
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} "mm10" &
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Reads_Profile_Promoter 'TSS' ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Peaks_Distribution_Analysis ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_Venn_Diagram ${__EXE_PATH} 'bed'
	#RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[i]} "mm10"
	#RUN_RPKM ${__INPUT_SAMPLE_DIR_List[i]} 'bed'
	#RUN_Reads_Profile_Promoter_genebody ${__INPUT_SAMPLE_DIR_List[i]}
	#RUN_TOPHAT ${__INPUT_SAMPLE_DIR_List[i]} "Treg" "mm10" "WT_Online_Ref"
	#RUN_BED2WIG ${__INPUT_SAMPLE_DIR_List[i]} ${SPECIES}
	
	#RUN_Wig2BigWig ${__RAW_DATA_PATH_DIR} ${__INPUT_SAMPLE_DIR_List[i]} 'Tcf1' ${SPECIES} ${Data_Provider}
	
	#RUN_CUFFDIFF ${__INPUT_SAMPLE_DIR_List[*]}
	break
#### FOR a full cycle, it must be clear its READS_DIR in the end.
	#echo "Unset DIR sets."
	#unset ${__FASTQ_DIR_R1} ${__FASTQ_DIR_R2}
	
	#### Parallel TEST
done
	

	echo "End Date: `date`"
	echo -e "\a FINISHED ALERT !"
	
	if [ ${Alert_email} == 0 ];then
	EMAIL_ME "${Start_Date}" ${Process_NAME}
	fi
}

##### Following Line is very IMPORTANT 
main "$@"

