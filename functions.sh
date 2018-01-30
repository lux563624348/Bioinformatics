#!/bin/bash
# /mnt/c/Users/Xiang/Dropbox/AResearch/version.beta

#set -e	#### terminating the script if any command exited with a nonzero exit status
set -u #### prevents the error by aborting the script if a variable’s value is unset
set -o pipefail 	#### check on the p398 of book Bioinformatics Data Skills.
#set -o noclobber ####The noclobber option tells bash not to overwrite any existing files when you redirect output.

## INTRODUCTION
########################################################################
## 12/01/2017
## By Xiang Li, basing on Zhouhao's scripts.
## lux@gwu.edu
## Peng's Lab
## Ver.1.0
########################################################################

## DECLARATION
########################################################################
  #  means NEED MODIFICATION. "VERY IMPORTANT INFORMATION"
  ## means Title_level 1.
 ### means Title_level 2.
#### means comment.
########################################################################

##	FUNDEMENTAL FUNCTIONS
########################################################################
RUN_TEST(){
	local AAA=$1
	local BBB=$2
	CHECK_arguments $# 2
	echo "function AAA is $AAA"
	read -p "Input a number from keyboard." Read_key
	#echo "Your input is: ${Read_key}"
	printf "Your input is: %s ${Read_key} \n"
} 

RUN_CHOOSE(){
	local default="[0 1 NO YES]"
	local answer=""
	
	read -p "Cancel Email Alert?" answer #$0 > $1 
	#### read -s (silent input)
	printf "%b" "\n"
	local iteration_num=0
	
	while [[ ${iteration_num} -lt 3 ]]
	do
		local refer=$(grep -io $answer <<< ${default}) ###-o only match part -i ignore-case
		if [ -z ${refer} ];then
		echo "Unexpected Answer ${answer}"
		read -p "Try Again. (Yes or No)  " answer
		iteration_num=$(expr ${iteration_num} + 1)
		elif [[ $(grep $refer <<< "1YES") = "1YES" ]];then
		echo "Email Alert Confirmed."
		return 1 
		break
		else
		echo "No Alert Confirmed."
		return 0
		break
		fi
	done
}

RUN_Download (){
### Download Login information and the download directory.
	web_address="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX184/SRX184860"
	user="xue"
	password="ayooxuxue"
	cd $__RAW_DATA_PATH_DIR
	DOWNLOAD_TARGET_NAME=$1;
########################################################################
	for ((i = 0; i <= $(expr $SAMPLE_NUM - 1); i++))
	do
####If download file is a folder. IT MUST END WITH trailing slash  "/"
	Down_dir=$__RAW_DATA_PATH_DIR/$DOWNLOAD_TARGET_NAME/;
	echo "wget --no-check-certificate -r -c -nH --user=$user --password=$pass_word --accept=gz --no-parent $web_address/${INPUT_SAMPLE_DIR_List[$i]}/"
	wget --no-check-certificate -r -c -nH --user=$user --password=$password --accept=gz --no-parent $web_adress
	echo "$DOWNLOAD_TARGET_NAME Downloaded Completed"
	echo ""
	done
	}

PRE_READS_DIR(){
	### PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	CHECK_arguments $# 2
	echo "Entering one Library of RAW_DATA_DIR: $__RAW_DATA_PATH_DIR/$1"
	echo "Searching file type of $2"
	local Current_DATA_DIR="${__RAW_DATA_PATH_DIR}/$1"
	cd ${Current_DATA_DIR}
	
	
########################################################################
	#local DIRLIST_R1="$( find -name "*R1.fastq" | sort -n | xargs)"
	#local DIRLIST_R2="$( find -name "*R2.fastq" | sort -n | xargs)" #
	
	local DIRLIST_R1_gz="$( find -name "*R1_goodreads.$2" | sort -n | xargs)"
	local DIRLIST_R2_gz="$( find -name "*R2_goodreads.$2" | sort -n | xargs)"


#### R1 Saving		
	local k=0
	for FILE_DIR in $DIRLIST_R1_gz
	do
		__FASTQ_DIR_R1[k]="${Current_DATA_DIR}/${FILE_DIR: 2}"
		echo "Saving R1 reads DIR as ${__FASTQ_DIR_R1[k]}"
		k=`expr $k + 1`
	done
	
#### R2 Saving	
	local k=0
	__FASTQ_DIR_R2[k]=""
	for FILE_DIR in $DIRLIST_R2_gz
	do
		__FASTQ_DIR_R2[k]="${Current_DATA_DIR}/${FILE_DIR: 2}"
		echo "Saving R2 reads DIR as ${__FASTQ_DIR_R2[k]}"
		k=`expr $k + 1`
	done
	
	echo "Finish Preparing READS of Library: $1"
	unset k
	}

RUN_SRA2FASTQ(){
	## Usage : RUN_SRA2FASTQ $sra_files_path
	echo "SRA to FASTQ"
	echo "fastq-dump -I --split-files --gzip $1"
	fastq-dump -I --split-files --gzip $1 
	}

RUN_READLINE(){
	#### USAGE: RUN_READLINE $INPUT $OUTPUT
	CHECK_arguments $# 2
	echo "READ_LINE"
	while read -a line
	do
		if [ ${line[6]} \> 4.0 ]; then
			#echo ${line[6]}
			echo ${line[*]} >> $2 # OUTPUT file
		fi
	done < $1   #INPUT file
	} 

RUN_CHANGE_STR(){
	#!/usr/bin/env bash
# cookbook filename: suffixer
#
# rename files that end in .bad to be .bash
for FN in *.bad
do
mv "${FN}" "${FN%bad}bash"
done
	}

##	MODOLES
########################################################################

RUN_FAST_QC_BT (){
####	RUN_FAST_QC ${fastq_file_path}
####	Usage: RUN_FAST_QC $INPUT_DATA_DIR
########################################################################
### FASTQC Output DIR setting ...
	local Output_fastqc=$__EXE_PATH/fastqc
	DIR_CHECK_CREATE ${Output_fastqc}
	cd $__RAW_DATA_PATH_DIR
	
	echo "fastqc -o ${Output_fastqc} $1"
	fastqc -o ${Output_fastqc} $1
	
	echo "$1 Fastqc Completed!"
	}

RUN_BOWTIE2(){
	### RUN_BOWTIE2 ${__INPUT_SAMPLE_DIR_List[3]}
	echo "RUN_BOWTIE2"
	CHECK_arguments $# 1
	BOWTIEINDEXS_mm9="/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
	
	#BOWTIEINDEXS_mm9_pMXs_combo="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/mm9_pMXs_combo/mm9_pMXs_combo"
	#BOWTIEINDEXS_mm9_combo_MMLV_pMXs="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/mm9_comdo_MMLV_pMXs/mm9_comdo_MMLV_pMXs"
	BOWTIEINDEXS_MMLV_pMXs_combo="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/MMLV_pMXs_combo/MMLV_pMXs_combo"
	BOWTIEINDEXS_MMLV="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/MMLV_index/MMLV"

	local INPUT_NAME=$1
	local BOWTIEINDEXS_name="BOWTIEINDEXS_mm9"
	local BOWTIEINDEXS="${BOWTIEINDEXS_mm9}"
	echo "Bowtie2 ref= $BOWTIEINDEXS_name"

	local SPECIES="mm9"
	local THREADS=4
	local OUTPUT_BOWTIE2_FOLDER="${__EXE_PATH}/${INPUT_NAME}/bowtie2_results"

	DIR_CHECK_CREATE ${OUTPUT_BOWTIE2_FOLDER}
########################################################################
	local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
	echo ""

	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	echo "Pair End Mode"
	
	 echo "bowtie2 -p $THREADS -t --no-unal --non-deterministic -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S $OUTPUT_BOWTIE2"
	 bowtie2 -p $THREADS -t --no-unal --non-deterministic -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S $OUTPUT_BOWTIE2
	
	#### concordantly pair output
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS_MMLV -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS_MMLV -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
	#### using un concordantly pair-ends do bowtie2 again.
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS_MMLV_pMXs_combo -1 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R1.fastq.gz -2 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R2.fastq.gz -S $OUTPUT_BOWTIE2"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS_MMLV_pMXs_combo -1 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R1.fastq.gz -2 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R2.fastq.gz -S $OUTPUT_BOWTIE2
	
	else
	echo "Single End Mode."
	echo "bowtie2 -p $THREADS --no-unal --non-deterministic -x $BOWTIEINDEXS -U ${__FASTQ_DIR_R1[0]} -S $OUTPUT_BOWTIE2"
	bowtie2 -p $THREADS --no-unal --non-deterministic -x $BOWTIEINDEXS -U ${__FASTQ_DIR_R1[0]} -S $OUTPUT_BOWTIE2
	fi

echo ""
###sam2bed
echo "samtools view $OUTPUT_BOWTIE2 -Sb | bamToBed -i stdin > $OUTPUT_BOWTIE2_FOLDER/$INPUT_NAME.bed"
samtools view $OUTPUT_BOWTIE2 -Sb | bamToBed -i stdin > $OUTPUT_BOWTIE2_FOLDER/$INPUT_NAME.bed 

echo ""
###bed2wigs
echo "RUN_bed2wig $OUTPUT_BOWTIE2_FOLDER $INPUT_NAME"
RUN_bed2wig $OUTPUT_BOWTIE2_FOLDER $INPUT_NAME

echo ""
###Wig2Bigwig
echo "RUN_Wig2BigWig $OUTPUT_BOWTIE2_FOLDER $INPUT_NAME"
RUN_Wig2BigWig $OUTPUT_BOWTIE2_FOLDER $INPUT_NAME

	}

RUN_bed_intersect(){
	#### usage RUN_bed_intersect $1 $2
	CHECK_arguments $# 2
	FILE_TYPE='bed'
	local Input_A1=${__RAW_DATA_PATH_DIR}/${1}.${FILE_TYPE}
	local Input_B1=${__RAW_DATA_PATH_DIR}/${2}.${FILE_TYPE}
	DIR_CHECK_CREATE "${__EXE_PATH}/Common"
	DIR_CHECK_CREATE "${__EXE_PATH}/Solo"
	
	local Title=""
	
	####################################################################
	echo "$Title"
	echo ""
	
	echo "Input A1"
	wc -l ${Input_A1}
	echo "Input B1"
	wc -l ${Input_B1}
	
	Output_Path="${__EXE_PATH}/Common/Common_${1}_vs_${2}"
	Output_Path2="${__EXE_PATH}/Solo/Solo_${1}_vs_${2}"

	echo "bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.bed"
	bedtools intersect -wa -u -a ${Input_A1} -b ${Input_B1} > ${Output_Path}.${FILE_TYPE}
	echo "bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.bed"
	bedtools intersect -v -a ${Input_A1} -b ${Input_B1} > ${Output_Path2}.${FILE_TYPE}
	
	echo "Output_Common "
	wc -l ${Output_Path}.${FILE_TYPE}
	echo "Output_Solo"
	wc -l ${Output_Path2}.${FILE_TYPE}
	
	}

RUN_bed2wig(){
	#### Usage: RUN_bed2wig <Output_DIR> <Output_name>
local SICER="/home/zzeng/Software/SICER1.1/SICER"
local EXEDIR="${SICER}/extra/tools/wiggle"
local WINDOW_SIZE=200
local FRAGMENT_SIZE=150
local OUTPUTDIR=$1
local OUTPUT=$2

sh $EXEDIR/bed2wig.sh $OUTPUTDIR $OUTPUT $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES
	}

RUN_SICER(){
#### Usage: RUN_SICER $1 $2 ($1 is input for SICER. $2 is the CONTRO Library)
echo "RUN_SICER"
CHECK_arguments $# 2
local EXEDIR="/home/lxiang/Software/SICER1.1/SICER"

#HIS_MODS_SET=(H3K27Ac H3K27Me3 H3K4Me1 H3K4Me3)
#GAP_SET=(400 600 400 200)

local GAP_SET=(400) # 600 400 200)
local REDUNDANCY=1
local INPUT_NAME=$1
local INPUI_CON=$2
local FRAGMENT_SIZE=150
local WINDOW_SIZE=200
local EFFECTIVEGENOME=0.85 
local FDR=0.0001





local IN_SICER_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
local IN_SICER_FILES=${INPUT_NAME}.bed

local CONTRO_SICER_DIR=${__EXE_PATH}/${INPUI_CON}/bowtie2_results
local CONTRO_SICER_FILE=${INPUI_CON}.bed

local OUT_SICER_FOLDER=${__EXE_PATH}/${INPUT_NAME}/SICER_results
DIR_CHECK_CREATE ${OUT_SICER_FOLDER}

cd ${OUT_SICER_FOLDER}

echo "sh ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} mm9 ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log"
sh ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} mm9 ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log



	}

RUN_MACS(){
#### Usage: RUN_MACS $1 $2 ($1 is input for MACS. $2 is the CONTRO Library)
echo "RUN_MACS"
CHECK_arguments $# 2
local EXEDIR="/home/lxiang/Software/python_tools/MACS-1.4.2/bin"

local INPUT_NAME=$1
local INPUI_CON=$2

local p_value=0.00001

local IN_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
local IN_FILES=${INPUT_NAME}.bam

local CONTRO_DIR=${__EXE_PATH}/${INPUI_CON}/bowtie2_results
local CONTRO_FILE=${INPUI_CON}.bam

local OUT_FOLDER=${__EXE_PATH}/${INPUT_NAME}/MACS_results
DIR_CHECK_CREATE ${OUT_FOLDER}

cd ${OUT_FOLDER}

echo "python ${EXEDIR}/macs14 --format BAMPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} -g 'mm' -n ${INPUT_NAME} -w -p ${p_value}"
python ${EXEDIR}/macs14 --format BAMPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} -g 'mm' -n ${INPUT_NAME} -w -p ${p_value} 
	}

RUN_MACS2(){
#### Usage: RUN_MACS2 $1 $2 ($1 is input for MACS. $2 is the CONTRO Library)
echo "RUN_MACS"
CHECK_arguments $# 2
local EXEDIR="/home/lxiang/Software/python_tools/MACS2-2.1.1.20160309/bin"

local INPUT_NAME=$1
local INPUI_CON=$2

#local FDR=0.05
local p_value=0.00001

local IN_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
local IN_FILES=${INPUT_NAME}.bed

local CONTRO_DIR=${__EXE_PATH}/${INPUI_CON}/bowtie2_results
local CONTRO_FILE=${INPUI_CON}.bed

local OUT_FOLDER=${__EXE_PATH}/${INPUT_NAME}/MACS2_results
DIR_CHECK_CREATE ${OUT_FOLDER}

echo "python ${EXEDIR}/macs2 callpeak --format BEDPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value}" # -m 4  -q ${FDR}"
python ${EXEDIR}/macs2 callpeak --format BEDPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value} # -m 4 -q ${FDR}

}

RUN_MACS2_Diff(){
#### Usage: RUN_MACS2_Diff $1 $2 $3 $4 $5 $6 $threshold ($1 is input for con1. $2 is the CONTRO_1  $3 is Con2. $4 is the CONTRO_2 )
echo "RUN_MACS2 Diff"
CHECK_arguments $# 7
local EXEDIR="/home/lxiang/Software/python_tools/MACS2-2.1.1.20160309/bin"

local CON1_NAME=$1
local CON1_TREAT=$2
local CON1_CONTRO=$3
local CON2_NAME=$4
local CON2_TREAT=$5
local CON2_CONTRO=$6

local p_value=0.00001
########################################################################
	if [ false ]; then
## Customized Part DIR
########################################################################
local CON1_FOLDER=${__EXE_PATH}/bowtie2_results
local CON1_FILE=${CON1_TREAT}.bed

local CON1_CONTRO_FOLDER=${__EXE_PATH}/bowtie2_results
local CON1_CONTRO_FILE=${CON1_CONTRO}.bed

local CON2_FOLDER=${__EXE_PATH}/bowtie2_results
local CON2_FILE=${CON2_TREAT}.bed

local CON2_CONTRO_FOLDER=${__EXE_PATH}/bowtie2_results
local CON2_CONTRO_FILE=${CON2_CONTRO}.bed
	fi

########################################################################
########################################################################
local CON1_FOLDER=${__EXE_PATH}/${CON1_TREAT}/bowtie2_results
local CON1_FILE=${CON1_TREAT}.bed

local CON1_CONTRO_FOLDER=${__EXE_PATH}/${CON1_CONTRO}/bowtie2_results
local CON1_CONTRO_FILE=${CON1_CONTRO}.bed

local CON2_FOLDER=${__EXE_PATH}/${CON2_TREAT}/bowtie2_results
local CON2_FILE=${CON2_TREAT}.bed

local CON2_CONTRO_FOLDER=${__EXE_PATH}/${CON2_CONTRO}/bowtie2_results
local CON2_CONTRO_FILE=${CON2_CONTRO}.bed

########################################################################

########################################################################
local OUT_FOLDER=${__EXE_PATH}/MACS2_Diff_results/${CON1_NAME}_vs_${CON2_NAME}
DIR_CHECK_CREATE ${OUT_FOLDER}

echo "python ${EXEDIR}/macs2 callpeak --format BEDPE -B -t ${CON1_FOLDER}/${CON1_FILE} -c ${CON1_CONTRO_FOLDER}/${CON1_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON1_NAME} -g 'mm' --nomodel --extsize 200"
python ${EXEDIR}/macs2 callpeak --format BEDPE -B -t ${CON1_FOLDER}/${CON1_FILE} -c ${CON1_CONTRO_FOLDER}/${CON1_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON1_NAME} -g 'mm' -p ${p_value} # --nomodel --extsize 200 

echo "python ${EXEDIR}/macs2 callpeak --format BEDPE -B -t ${CON2_FOLDER}/${CON2_FILE} -c ${CON2_CONTRO_FOLDER}/${CON2_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON2_NAME} -g 'mm' --nomodel --extsize 200"
python ${EXEDIR}/macs2 callpeak --format BEDPE -B -t ${CON2_FOLDER}/${CON2_FILE} -c ${CON2_CONTRO_FOLDER}/${CON2_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON2_NAME} -g 'mm' -p ${p_value} # --nomodel --extsize 200 

cd ${OUT_FOLDER}

local Size_Treat1=$(egrep "fragments after filtering in treatment" ${CON1_NAME}_peaks.xls | grep -o '[0-9]*')
local Size_Contro1=$(egrep "fragments after filtering in control" ${CON1_NAME}_peaks.xls | grep -o '[0-9]*')

local Size_Treat2=$(egrep "fragments after filtering in treatment" ${CON2_NAME}_peaks.xls | grep -o '[0-9]*')
local Size_Contro2=$(egrep "fragments after filtering in control" ${CON2_NAME}_peaks.xls | grep -o '[0-9]*')

if [ ${Size_Treat1} -lt ${Size_Contro1} ];then 
	local Size_1=${Size_Treat1}
	else
	local Size_1=${Size_Contro1}
	fi
	
if [ ${Size_Treat2} -lt ${Size_Contro2} ];then 
	local Size_2=${Size_Treat2}
	else
	local Size_2=${Size_Contro2}
	fi
	
	echo "Size_Treat1: ${Size_Treat1} and Size_Contro1: ${Size_Contro1}, Size_1 equal ${Size_1}"
	echo "Size_Treat2: ${Size_Treat2} and Size_Contro2: ${Size_Contro2}, Size_2 equal ${Size_2}"

macs2 bdgdiff --t1 ${CON1_NAME}_treat_pileup.bdg --c1 ${CON1_NAME}_control_lambda.bdg --t2 ${CON2_NAME}_treat_pileup.bdg \
--c2 ${CON2_NAME}_control_lambda.bdg --d1 ${Size_1} --d2 ${Size_2} -g 60 -l 120 --o-prefix diff_${CON1_NAME}_vs_${CON2_NAME}

}

RUN_CHECK_SEQ(){
#### Check a specific seq in fastq.gz
####	1.Raw_data 2.Seq
	cd $1   #### Entering the Raw_Data_DIR 
	OUT_REPORT="$1/seq_check_report"
	DIR_CHECK_CREATE $OUT_REPORT
	CHECK_SEQ=$2
	
	filename="${OUT_REPORT}/same_seq_counts.txt"
#### Head for oupput report	
	echo "Checing Sequence: $CHECK_SEQ " >>$filename
	echo "SAME_Seq_LINES  TOTAL_LINES  FILE_NAME" >>$filename
	
	echo "Detecting all files ..."
	DIRLIST_R="$( find -name "*.fastq.gz"| xargs)"
	echo "Saving the DIR of all fastq.gz files ..."

	for FILE_DIR in $DIRLIST_R
	do
	echo "gunzip -c ${FILE_DIR: 2}"
	gunzip ${FILE_DIR: 2}
####	Check GATCGGAAGCG
	TOTAL_LINES=$(wc -l ${FILE_DIR: 2: -3})
	LINES=$(grep $CHECK_SEQ ${FILE_DIR: 2: -3} | wc -l)
	
	echo ""
	echo "$LINES	$TOTAL_LINES" >>$filename
	echo "gzip ${FILE_DIR: 2: -3}"
	
	gzip ${FILE_DIR: 2: -3}
	
	done
}

RUN_Counting_READS_Single_End(){
	local INPUT_NAME=$1
	cd $__EXE_PATH/$INPUT_NAME/bowtie2_results
	
	local NUM_Klf4=$( grep Klf4 $INPUT_NAME.sam | wc -l )
	local NUM_Oct4=$( grep Oct4 $INPUT_NAME.sam | wc -l )
	local NUM_Sox2=$( grep Sox2 $INPUT_NAME.sam | wc -l )
	local NUM_chr=$( grep chr $INPUT_NAME.sam | wc -l )
	
	local chr_multiple_Klf4=$( grep chr $INPUT_NAME.sam | grep Klf4 | wc -l )
	local Oct4_multiple_Klf4=$( grep Oct4 $INPUT_NAME.sam | grep Klf4 | wc -l )
	local Sox2_multiple_Klf4=$( grep Sox2 $INPUT_NAME.sam | grep Klf4 | wc -l )
	
	local chr_multiple_Oct4=$( grep chr $INPUT_NAME.sam | grep Oct4 | wc -l )
	local Klf4_multiple_Oct4=$( grep Klf4 $INPUT_NAME.sam | grep Oct4 | wc -l )
	local Sox2_multiple_Oct4=$( grep Sox2 $INPUT_NAME.sam | grep Oct4 | wc -l )
	
	
	local chr_multiple_Sox2=$( grep chr $INPUT_NAME.sam | grep Sox2 | wc -l )
	local Oct4_multiple_Sox2=$( grep Oct4 $INPUT_NAME.sam | grep Sox2 | wc -l )
	local Klf4_multiple_Sox2=$( grep Klf4 $INPUT_NAME.sam | grep Sox2 | wc -l )
	
	
	local NUM_chr=$(expr $NUM_chr - $chr_multiple_Klf4 - $chr_multiple_Oct4 - $chr_multiple_Sox2 - 22 )
	local NUM_Klf4=$(expr $NUM_Klf4 - $Oct4_multiple_Klf4 - $Sox2_multiple_Klf4 - $chr_multiple_Klf4 - 1 )
	local NUM_Oct4=$(expr $NUM_Oct4 - $chr_multiple_Oct4 - $Klf4_multiple_Oct4 - $Sox2_multiple_Oct4 - 1 )
	local NUM_Sox2=$(expr $NUM_Sox2 - $chr_multiple_Sox2 - $Oct4_multiple_Sox2 - $Klf4_multiple_Sox2 - 1 )
	
	local filename=$__EXE_PATH/counting_results.log
	
	if [ ! -f $filename ];then
	echo "The counting results of uniquely mapping is:" >$filename
	echo "Counting file: $INPUT_NAME">>$filename
	echo "chr: $NUM_chr" >>$filename
	echo "pMXs-Klf4: $NUM_Klf4" >>$filename
	echo "pMXs-Oct4: $NUM_Oct4" >>$filename
	echo "pMXs-Sox2: $NUM_Sox2" >>$filename
	else
	echo ""
	echo "Counting file: $INPUT_NAME">>$filename
	echo "chr: $NUM_chr" >>$filename
	echo "pMXs-Klf4: $NUM_Klf4" >>$filename
	echo "pMXs-Oct4: $NUM_Oct4" >>$filename
	echo "pMXs-Sox2: $NUM_Sox2" >>$filename
	fi
	
	}

RUN_Counting_READS_Pair_End(){
	local INPUT_NAME=$1
	cd $__EXE_PATH/$INPUT_NAME/bowtie2_results
	
	local Q_NAME=HS2
	#local Q_NAME=HWI
	
	local NUM_Klf4=$( grep ${Q_NAME} pMXs-Klf4-full_R1.fastq | wc -l )
	local NUM_Oct4=$( grep ${Q_NAME} pMXs-Oct4-full_R1.fastq | wc -l )
	local NUM_Sox2=$( grep ${Q_NAME} pMXs-Sox2-full_R1.fastq | wc -l )
	local NUM_MMLV=$( grep ${Q_NAME} MMLV_R1.fastq | wc -l )
	
	local NUM_MMLV=$(expr $NUM_MMLV \* 2)
	local NUM_Klf4=$(expr $NUM_Klf4 \* 2)
	local NUM_Oct4=$(expr $NUM_Oct4 \* 2)
	local NUM_Sox2=$(expr $NUM_Sox2 \* 2)
	
	local filename=$__EXE_PATH/counting_results.log
	
	if [ ! -f $filename ];then
	echo "The counting results of uniquely mapping is:" >$filename
	echo "Counting file: $INPUT_NAME">>$filename
	echo "MMLV: $NUM_MMLV" >>$filename
	echo "pMXs-Klf4: $NUM_Klf4" >>$filename
	echo "pMXs-Oct4: $NUM_Oct4" >>$filename
	echo "pMXs-Sox2: $NUM_Sox2" >>$filename
	else
	echo ""
	echo "Counting file: $INPUT_NAME">>$filename
	echo "MMLV: $NUM_MMLV" >>$filename
	echo "pMXs-Klf4: $NUM_Klf4" >>$filename
	echo "pMXs-Oct4: $NUM_Oct4" >>$filename
	echo "pMXs-Sox2: $NUM_Sox2" >>$filename
	fi
	
	if [ -e *.fastq ];then
	gzip *.fastq
	fi
	}

RUN_CUFFDIFF(){
	
OUTPUT_tophat=$1/tophat_anlysis
OUT_SAMPLE_NAME=$2
DIR_CHECK_CREATE $OUTPUT_tophat/cuffdiff
SAMPLE_NUM=${#OUT_SAMPLE_NAME{*}}

#### Preparing the .bam files for cuffdiff
for ((i = 0; i <= $(expr $SAMPLE_NUM - 1); i++))
do	
	DATA_BAM[$i]=$OUTPUT_tophat/${OUT_SAMPLE_NAME[i]}/${OUT_SAMPLE_NAME[i]}.bam
	echo ${DATA_BAM[$i]}
done

### TWO GROUP, Four DATA FILES, 
########################################################################
if [ $SAMPLE_NUM -eq "4" ]
then
	echo "$SAMPLE_NUM Data files are loading..."
	echo "cuffdiff -o $OUTPUTDIR_CuffDiff -p $THREADS -L ${SAMPLENAME_GROUP[0]},${SAMPLENAME_GROUP[1]} $GTFFILE ${DATA_BAM[0]},${DATA_BAM[1]} ${DATA_BAM[2]},${DATA_BAM[3]}"
	cuffdiff -o $OUTPUTDIR_CuffDiff -p $THREADS -L ${SAMPLENAME_GROUP[0]},${SAMPLENAME_GROUP[1]} $GTFFILE ${DATA_BAM[0]},${DATA_BAM[1]} ${DATA_BAM[2]},${DATA_BAM[3]}
fi
########################################################################

echo "CuffDiff Completed!"
	}

RUN_FAST_QC (){
####	RUN_FAST_QC ${fastq_file_path}
####	Usage: RUN_FAST_QC $INPUT_DATA_DIR
########################################################################
### FASTQC Output DIR setting ...
	local Output_fastqc=$__EXE_PATH/fastqc
	DIR_CHECK_CREATE ${Output_fastqc}
	cd $__RAW_DATA_PATH_DIR
	
	for fastq_file in $__FASTQ_DIR_R1
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} $fastq_file
	done
	
	for fastq_file in $__FASTQ_DIR_R2
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} $fastq_file
	done
	
	echo "Fastq Completed!"
	}

RUN_REMOVE_SEQ(){
## USAGE:	RUN_REMOVE_SEQ ${__INPUT_SAMPLE_DIR_List[i]}
#### Check a specific seq in fastq.gz
####	1.Raw_data 2.Seq
local bbmap=/home/lxiang/Software/bbmap
local REMOVE_SEQ=/home/lxiang/Raw_Data/Paul/34bc/Reference_gene_seq/MLVs/MMLV.fa
local OUT=$__EXE_PATH/$1
DIR_CHECK_CREATE $OUT
local KERS=$2
########################################################################

echo "nohup $bbmap/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stas.txt"
nohup $bbmap/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stats.txt

########################################################################
	}

RUN_SELECT_SEQ_SAM(){
	#### Usage: RUN_SELECT_SEQ_SAM $INPUT_SAM_FOLDER $INPUT_NAME $RNAME
	#### Given a sam file DIR
	local INPUT_NAME=$1
	local INPUT_SAM_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
	local OUTPUT_RNAME=$2
	CHECK_arguments $# 2
	
	echo "Entering one Library of SAM_DATA_FOLDER: $INPUT_SAM_FOLDER"
	cd ${INPUT_SAM_FOLDER}
	local INPUT_SAM=${INPUT_NAME}.sam
	local SORT_BAM=${INPUT_NAME}\_sorted.bam
	#local SORT_SAM=$INPUT_NAME\_sorted.sam
	
	if [ ! -e $SORT_BAM ];then
	echo "samtools sort -@ 4 -o ${SORT_BAM} -O BAM ${INPUT_SAM}"
	samtools sort -@ 4 -o ${SORT_BAM} -O BAM ${INPUT_SAM}
	
	echo "Create a BAI INDEX: samtools index -b ${SORT_BAM}"
	samtools index -b ${SORT_BAM}
	fi 
	
	echo "Only output alignment with the ref: ${OUTPUT_RNAME}"
	samtools view -h -b ${SORT_BAM} ${OUTPUT_RNAME} -o ${OUTPUT_RNAME}.bam
	
	echo "samtools view -@ 4 -h ${SORT_BAM} ${OUTPUT_RNAME} -o ${OUTPUT_RNAME}.sam"
	samtools view -h ${SORT_BAM} ${OUTPUT_RNAME} -o ${OUTPUT_RNAME}.sam
	
	echo "samtools sort -n -o ${OUTPUT_RNAME}.bam -O BAM ${OUTPUT_RNAME}.sam"
	samtools sort -@ 4 -n -o ${OUTPUT_RNAME}.bam -O BAM ${OUTPUT_RNAME}.sam
	samtools sort -@ 4 -n -o ${OUTPUT_RNAME}.sam -O SAM ${OUTPUT_RNAME}.sam
	
	#### Single End reads
	#bedtools bamtofastq -i $OUTPUT_RNAME.bam -fq $INPUT_SAM_FOLDER/$OUTPUT_RNAME\_R1.fastq
	
	#### Pair Ends reads
	echo "bedtools bamtofastq -i ${OUTPUT_RNAME}.bam  -fq ${INPUT_SAM_FOLDER}/${OUTPUT_RNAME}\_R1.fastq  -fq2 ${INPUT_SAM_FOLDER}/${OUTPUT_RNAME}\_R2.fastq"
	bedtools bamtofastq -i ${OUTPUT_RNAME}.bam  -fq ${INPUT_SAM_FOLDER}/${OUTPUT_RNAME}\_R1.fastq  -fq2 ${INPUT_SAM_FOLDER}/${OUTPUT_RNAME}\_R2.fastq
	
	#rm *.bam

	}

RUN_GZIP(){
####	run gzip for all files under input DIR
####	Usage: 
########################################################################
### gzip
	cd ${__RAW_DATA_PATH_DIR}
########################################################################	
    echo "Convert all matched*.fastq to fastq.gz"
	local DIRLIST_R="$( find -name "*.fastq"| xargs)"

	for FILE_DIR in ${DIRLIST_R}
	do
	echo "gzip ${__RAW_DATA_PATH_DIR}/${FILE_DIR: 2}"
	gzip ${__RAW_DATA_PATH_DIR}/${FILE_DIR: 2}
	done
	}

RUN_Wig2BigWig (){
	###RUN_Wig2BigWig   $INPUT_FOLDER $INPUT_NAME
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $Sample_Wig_NAME
	CHECK_arguments $# 2
	local Data_provider=Haihui
	local Data_label=Lef1
		
	local Tracks_NAME=$2
	#local Contro_Tracks_NAME=$2
	local Tracks_NAME_Lable=${Tracks_NAME: 7:-12} ##Skip out of Sample_ (7) and _20160827000 (-12)
	#local Contro_Tracks_NAME_Lable=${Contro_Tracks_NAME: 7:-12}
	
	####INPUT
	local UCSC_DIR=/home/lxiang/Software/UCSC
	local Reference_genome_sizes=$UCSC_DIR/genome_sizes/mm9.chrom.sizes
	local Wig_DIR=$1
	cd ${Wig_DIR}
#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigWigs

#### Modified Tracks_NAME
	local Tracks_NAME=${Tracks_NAME}-W200-RPKM
	echo "$UCSC_DIR/wigToBigWig ${Tracks_NAME}.wig $Reference_genome_sizes ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	$UCSC_DIR/wigToBigWig ${Tracks_NAME}.wig $Reference_genome_sizes ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw

###	Creating hub.txt
	cd $OUTPUTDIR_tracks_hub
	local filename=hub.txt
	if [ ! -f $filename ];then
	echo "hub ${Data_label}" >>$filename
	echo "shortLabel ${Data_label}" >>$filename
	echo "longLabel ${Data_label}_/${Data_provider}" >>$filename
	echo "genomesFile genomes.txt" >>$filename
	echo "email lux@gwu.edu" >>$filename
	echo "descriptionUrl descriptionUrl" >>$filename
	fi
########################################################################

###	Creating genomes.txt	
	trackDb_NAME="mm9"
	local filename=genomes.txt
	if [ ! -f $filename ];then
	echo "genome $trackDb_NAME" >>$filename
	echo "trackDb $trackDb_NAME/trackDb.txt" >>$filename
	fi
########################################################################

	DIR_CHECK_CREATE $OUTPUTDIR_tracks_hub/${trackDb_NAME}	
	cd $OUTPUTDIR_tracks_hub/${trackDb_NAME}
	
	local filename=trackDb.txt
	local bw_Url=https://s3.amazonaws.com/xianglilab/tracks_hub/${Data_provider}/${Data_label}/BigWigs/${Tracks_NAME}.bw
	echo "track ${Tracks_NAME}" >>$filename
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "visibility full" >>$filename
	echo "maxHeightPixels 60" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
	}

RUN_BigGraph2BigWig (){
	##RUN_BigGraph2BigWig   $INPUT_FOLDER $INPUT_NAME
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $Sample_Wig_NAME
	CHECK_arguments $# 2
	local Data_provider=Haihui
	local Data_label=Lef1_MACS2
		
	local Tracks_NAME=$2
	#local Contro_Tracks_NAME=$2
	local Tracks_NAME_Lable=${Tracks_NAME: 7:-12} ##Skip out of Sample_ (7) and _20160827000 (-12)
	#local Contro_Tracks_NAME_Lable=${Contro_Tracks_NAME: 7:-12}
	
	####INPUT
	local UCSC_DIR=/home/lxiang/Software/UCSC
	local Reference_genome_sizes=$UCSC_DIR/genome_sizes/mm9.chrom.sizes
	local bdg_DIR=$1
	cd ${bdg_DIR}
#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigWigs
	
########################################################################
#### Modified Tracks_NAME
	#local Tracks_NAME=${Tracks_NAME}_treat_pileup
	##control sample
	local Tracks_NAME=${Tracks_NAME}_control_lambda
	
	
	if [ ! -f ${Tracks_NAME}_sorted.bdg ];then
	echo "sort -k1,1 -k2,2n ${Tracks_NAME}.bdg > ${Tracks_NAME}_sorted.bdg"
	sort -k1,1 -k2,2n ${Tracks_NAME}.bdg > ${Tracks_NAME}_sorted.bdg
	fi
	
	echo "$UCSC_DIR/bedGraphToBigWig ${Tracks_NAME}_sorted.bdg $Reference_genome_sizes ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	$UCSC_DIR/bedGraphToBigWig ${Tracks_NAME}_sorted.bdg $Reference_genome_sizes ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw

###	Creating hub.txt
	cd $OUTPUTDIR_tracks_hub
	local filename=hub.txt
	if [ ! -f $filename ];then
	echo "hub ${Data_label}" >>$filename
	echo "shortLabel ${Data_label}" >>$filename
	echo "longLabel ${Data_label}_/${Data_provider}" >>$filename
	echo "genomesFile genomes.txt" >>$filename
	echo "email lux@gwu.edu" >>$filename
	echo "descriptionUrl descriptionUrl" >>$filename
	fi
########################################################################

###	Creating genomes.txt	
	trackDb_NAME="mm9"
	local filename=genomes.txt
	if [ ! -f $filename ];then
	echo "genome $trackDb_NAME" >>$filename
	echo "trackDb $trackDb_NAME/trackDb.txt" >>$filename
	fi
########################################################################

	DIR_CHECK_CREATE $OUTPUTDIR_tracks_hub/${trackDb_NAME}	
	cd $OUTPUTDIR_tracks_hub/${trackDb_NAME}
	
	local filename=trackDb.txt
	local bw_Url=https://s3.amazonaws.com/xianglilab/tracks_hub/${Data_provider}/${Data_label}/BigWigs/${Tracks_NAME}.bw
	echo "track ${Tracks_NAME}" >>$filename
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "visibility full" >>$filename
	echo "maxHeightPixels 60" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
}

RUN_TOPHAT (){
#### Usage: RUN_TOPHAT 1.$RAW_DATA_PATH_DIR 2.$EXE_PATH 3.${INPUT_SAMPLE_DIR_List[*]} 4.OUT_SAMPLE_NAME_List
### Tophat Input and Output DIR setting ...	
SICER_DIR=/home/lxiang/Software/SICER1.1/SICER
EXEDIR=$SICER_DIR/extra/tools/wiggle
GTFFILE=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2014-05-23-16-05-24/Genes/genes.gtf
BOWTIEINDEXS=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
### Operation PARAMETERS Setting
THREADS=4
SPECIES=mm9
WINDOW_SIZE=10
FRAGMENT_SIZE=0
########################################################################
RAW_DATA_PATH_DIR=$1/$3
EXE_PATH_tophat_analysis=$2/tophat_analysis
INPUT_SAMPLE_DIR=$3
OUT_SAMPLE_DIR=$4
	
DIR_CHECK_CREATE $EXE_PATH_tophat_analysis
	
	echo "Entering RAW_DATA_DIR: $RAW_DATA_PATH_DIR"
	cd $RAW_DATA_PATH_DIR
	echo "Entered folder: $INPUT_SAMPLE_DIR"
	DIRLIST_R1="$( find -name "*R1*.gz" | xargs)"
	DIRLIST_R2="$( find -name "*R2*.gz" | xargs)"
#### R1 Saving		
	i=0
	for FILE_DIR in $DIRLIST_R1
	do
		FASTQ_DIR_R1[$i]=$RAW_DATA_PATH_DIR/${FILE_DIR: 2}
		echo "Saving R1 reads of ${FILE_DIR: 2}"
		i=`expr $i + 1`
	done
	
#### R2 Saving	
	i=0
	for FILE_DIR in $DIRLIST_R2
	do
		FASTQ_DIR_R2[$i]=$RAW_DATA_PATH_DIR/${FILE_DIR: 2}
		echo "Saving R2 reads of ${FILE_DIR: 2}"
		i=`expr $i + 1`
	done
	echo ""

#### Decide single end or pair ends mode
#### NOW it is only compatible with single file. Not with pieces files.
	if [ -n "${FASTQ_DIR_R1[0]}" -a -n "${FASTQ_DIR_R2[0]}" ]
	then
	echo "Pair End Mode"
	echo "tophat -p $THREADS -o $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR --GTF $GTFFILE $BOWTIEINDEXS ${FASTQ_DIR_R1[0]} ${FASTQ_DIR_R2[0]}"
	tophat -p $THREADS -o $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR --GTF $GTFFILE $BOWTIEINDEXS ${FASTQ_DIR_R1[0]} ${FASTQ_DIR_R2[0]}
	else
	echo "Single End Mode."
	echo "tophat -p $THREADS -o $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR --GTF $GTFFILE $BOWTIEINDEXS ${FASTQ_DIR_R1[0]}"
	tophat -p $THREADS -o $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR --GTF $GTFFILE $BOWTIEINDEXS ${FASTQ_DIR_R1[0]}
	fi

	echo "cd $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR"	
	cd $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR
	mv accepted_hits.bam $OUT_SAMPLE_DIR.bam
	
	echo "sh $EXEDIR/bam2wig.sh $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR $OUT_SAMPLE_DIR $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES"
	echo ""
	sh $EXEDIR/bam2wig.sh $EXE_PATH_tophat_analysis/$OUT_SAMPLE_DIR $OUT_SAMPLE_DIR $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES
	
	RUN_Wig2BigWig $OUT_SAMPLE_DIR
	
	rm $OUT_SAMPLE_DIR.bed
	echo "One tophat is completed."
	echo ""
	}

RUN_Peaks_Distribution_Analysis(){
	CHECK_arguments $# 1
	
	EXEDIR=/home/lxiang/Software/python_tools

	GTFDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gtf_files
	GTFFILE=mm9_genes.gtf

	PROMOTER_UPSTREAM_EXTENSION=1001
	TSS_REGION_LENGTH=2000


	INPUTFILE=${1}.bed
	OUTPUTFILE=${1}_distribution.txt

	echo "python $EXEDIR/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE"
	python $EXEDIR/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE
	echo ""
	echo ""

	}

RUN_CUT_Columns (){
	####RUN_CUT_Columns $Filename $start $end
	CHECK_arguments $# 3
	local File_type='bed'
	local File_Path=${__RAW_DATA_PATH_DIR}/${1}.${File_type}
	local start=${2}
	local end=${3}
	local delimiter="	"
	cut -f ${start}-${end} ${File_Path} > ${__RAW_DATA_PATH_DIR}/${1}_13.${File_type}
	echo ""
	echo ""
	}

########################################################################

##BASIC FUNCTIONS
########################################################################
EMAIL_ME(){
	CHECK_arguments $# 1
	echo "Start at $1 " | mailx -v -s "Project Finished" lux@gwu.edu
	}

CHECK_arguments(){
#### Usage: CHECK_arguments $# <NUM_arguments>
#### $# gives the number of command-line arguments.
#### <NUM_arguments> is the right number of arguments.
	if [ "$1" -lt $2 ];then 
	echo "error: too few arguments, you provided $1 arguments, $2 is required."
	exit 1
	fi
	} 

DIR_CHECK_CREATE (){
### Saving DIR Check and Create
	echo "Dir check and create are: $@"
	if [ ! -d $@ ];then
	mkdir -p $@
	fi
	}

BACKUP_FOLDER (){
### Saving DIR Check and Create
	echo "Dir Under backup are: $@"
	local $BP_DIR=back_up
	DIR_CHECK_CREATE ${BP_DIR}
	
	for filename in $@
	do
		cp ${filename} ${BP_DIR}/${filename}_bp
	done
	
	echo "Backup Finished."
	}


CLEAN(){
####Usage: RUN_CLEAN $Directory $File type.
	CHECK_arguments $# 2
	local Dir=$1
	local File_type=$2
	if [[ -n ${Dir} && -n ${File_type} ]]; then
	cd ${Dir}
		if (( ! $? )); then rm *${File_type}; fi
	else echo "ERR: Empty. Exit"
	fi
	}
########################################################################
Max(){
## Usage: Larger_value = $(RUN_Max $num_a $num_b)
	CHECK_arguments $# 2
	if [ $1 -gt $2 ]
	then
		echo $1
	else
		echo $2
	fi
	}


##### Following Line is very IMPORTANT 
#main "$@"
