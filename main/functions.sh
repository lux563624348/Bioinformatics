#!/bin/bash
# /mnt/c/Users/Xiang/Dropbox/AResearch/version.beta

#set -e	#### terminating the script if any command exited with a nonzero exit status
set +H #### to disable the C shell-style command history,‘!’ as a special character.
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
########################################################################
########################################################################

########################################################################
###    Global Variables
	#Shang
	UCSC_DIR=/opt/tools/UCSC
	Tools_DIR=/opt/tools
	Python_Tools_DIR=~/cloud_research/PengGroup/XLi/Python_tools
	Annotation_DIR=/home/Data/Annotation
	THREADS=16
	# Tang
	#Annotation_DIR=/home/Data/Annotation
	#THREADS=16
	#


######################################

##	FUNDEMENTAL FUNCTIONS FOR ANALYSIS MODULES
PRE_READS_DIR(){
	### PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz" "R1/R2 or Pairs"
	CHECK_arguments $# 3 # No less than 3
	echo "Entering one Library of RAW_DATA_DIR: $__RAW_DATA_PATH_DIR"
	echo "Searching file type as: $1 + $3 + $2"
	local Current_DATA_DIR=${__RAW_DATA_PATH_DIR}
	cd ${__RAW_DATA_PATH_DIR}
	Pair_Type=${3}
	
########################################################################
	#local DIRLIST_R1="$( find -name "*R1.fastq" | sort -n | xargs)"
	#local DIRLIST_R2="$( find -name "*R2.fastq" | sort -n | xargs)" #
	
	case ${Pair_Type} in
	"R1")
	echo "Single End, Only ${Pair_Type} "
	local DIRLIST_R1_gz="$( find -name "*${1}*R1*.${2}" | sort -n | xargs)"
	local DIRLIST_R2_gz=" "
	;;
	"R2")
	echo "Single End, Only ${Pair_Type} "
	local DIRLIST_R1_gz="$( find -name "*${1}*R2*.${2}" | sort -n | xargs)"
	local DIRLIST_R2_gz=" "
	;;
	"SRA")
	echo "SRA Mode"
	local DIRLIST_R1_gz="$( find -name "*${1}*_1*.${2}" | sort -n | xargs)"
	local DIRLIST_R2_gz="$( find -name "*${1}*_2*.${2}" | sort -n | xargs)"
	;;
	"Pairs")
	echo "Pair End, Both"
	local DIRLIST_R1_gz="$( find -name "*${1}*R1*.${2}" | sort -n | xargs)"
	local DIRLIST_R2_gz="$( find -name "*${1}*R2*.${2}" | sort -n | xargs)"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
	
	## for 34bc
	#local DIRLIST_R1_gz="$( find -name "un_conc_aligned_R1.$2" | sort -n | xargs)"
	#local DIRLIST_R2_gz="$( find -name "un_conc_aligned_R2.$2" | sort -n | xargs)"


#### R1 Saving		
	local k=0
	for FILE_DIR in ${DIRLIST_R1_gz[*]}
	do
		__FASTQ_DIR_R1[k]="${Current_DATA_DIR}/${FILE_DIR: 2}"
		echo "Saving R1 reads DIR as ${__FASTQ_DIR_R1[k]}"
		k=`expr $k + 1`
	done
	
#### R2 Saving	
	local k=0
	__FASTQ_DIR_R2[k]=""
	for FILE_DIR in ${DIRLIST_R2_gz[*]}
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
	#https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	CHECK_arguments $# 1
	## Given a path contains SraAccList.txt
	## Download all SRA files.
	local exe_path=/opt/tools/sratoolkit.2.9.2-ubuntu64/bin
	cd ${__RAW_DATA_PATH_DIR}
	#${exe_path}/prefetch --option-file SraAccList.txt
	##
	local INPUT_SRA_LIST=${1} #(cat SraAccList.txt)
	echo "SRA to FASTQ"
	for SRA in ${INPUT_SRA_LIST[*]}
	do
		echo "fastq-dump –X 5 –Z -split-files --gzip ${SRA}"
		${exe_path}/fastq-dump –X 5 –Z -split-files --gzip ${SRA}
	done
	}

RUN_COPY_AND_CHANGE_NAME(){
	### RUN_COPY_AND_CHANGE_NAME $1 $2 $3 ($1 is INPUT_NAME $2> OUTPUT_NAME $3 "fastq.gz")
	CHECK_arguments $# 3
	local INPUT_DATA_PATH="${__RAW_DATA_PATH_DIR}/${1}"
	local OUTPUT_DATA_PATH="${__EXE_PATH}/${2}"
	
	mkdir -p ${OUTPUT_DATA_PATH}
	
	cp ${INPUT_DATA_PATH}.${3} ${OUTPUT_DATA_PATH}/${2}.${3}
	echo "Finished Copying File: ${2}"
	
	
	
	}

RUN_BBDUK_Trimming (){
####	RUN_BBDUK_Trimming
####	Usage: RUN_BBDUK_Trimming $INPUT_DATA_DIR
########################################################################
### FASTQC Output DIR setting ...
	local Output_Trimmed=${__RAW_DATA_PATH_DIR}/Trimmed_Results
	DIR_CHECK_CREATE ${Output_Trimmed}
	cd ${__RAW_DATA_PATH_DIR}
	
	#Force-Trimming:
	#bbduk.sh in=reads.fq out=clean.fq ftl=10 ftr=139
	#  This will trim the leftmost 10 bases (ftl=10) and also trim the right end after to position 139 (zero-based). The resulting read would be 130bp long. For example, a 150bp read would have the first 10 bases trimmed (bases 0-9, keeping 10+) and the last 10 bases trimmed (bases 140-149, keeping 139 and lower).

	for file in ${__FASTQ_DIR_R1[*]}
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} ${fastq_file}
	done
	
	for fastq_file in ${__FASTQ_DIR_R2[*]}
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} ${fastq_file}
	done
	
	echo "Fastq Completed!"
	}

RUN_FAST_QC (){
####	RUN_FAST_QC
####	Usage: RUN_FAST_QC $INPUT_DATA_DIR
########################################################################
### FASTQC Output DIR setting ...
	local Output_fastqc=${__RAW_DATA_PATH_DIR}/fastqc
	DIR_CHECK_CREATE ${Output_fastqc}
	cd ${__RAW_DATA_PATH_DIR}
	
	for fastq_file in ${__FASTQ_DIR_R1[*]}
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} ${fastq_file} &
	done
	
	for fastq_file in ${__FASTQ_DIR_R2[*]}
	do
	echo "fastqc -o ${Output_fastqc} $fastq_file"
	fastqc -o ${Output_fastqc} ${fastq_file} &
	done
	
	echo "cd ${Output_fastqc}"
	cd ${Output_fastqc}
	
	#if [ -f multiqc_report.html ];then
	#echo "multiqc *fastqc.zip --ignore *.html"
	#multiqc *fastqc.zip --ignore *.html
	#fi

	
	echo "Fastq Completed!"
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

RUN_Read_Process_Data(){
	#### USAGE: RUN_READLINE $INPUT $OUTPUT
	CHECK_arguments $# 2
	echo "READ_LINE"
	while read -a line
	do
		if [ ${line[6]} \> 0.0 ]; then
			echo ${line[0]} ${line[1]} $(expr ${line[6]} + ${line[1]}) ${line[6]} >> $2 # OUTPUT file
		fi
	done < $1 #INPUT file
	}

RUN_awk_Process_Data(){
	#### USAGE: RUN_READLINE $INPUT $OUTPUT
	CHECK_arguments $# 2
	echo "READ_LINE"
	awk 'BEGIN{print "Don\47t Panic!"}'
		
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

RUN_Venn_Diagram(){
	#### Usage: RUN_Bedtools_Merge $1 $2
	### For bed format Venn
	:'For example:
	Chr	 start	 end
	File A: Chr1 	1	20
	File B: Chr1	15	30
	File C: Chr1	45	60
	> Union> File: union.bed
	Chr	 start	 end
	Chr1 	1	30
	Chr1	45	60

	Overlap union.bed with A, B and C one by one.
				‘A’	‘B’	‘C’
	Chr1 	1	30	1	1	0  
	Chr1	45	60	0	0	0

	Because chr1	1	30 has overlap both in ‘A’ and ‘B’, then we says that this peak is common in A and B.'
	
	CHECK_arguments $# 2
	local Current_DATA_DIR=${1}
	local FILE_TYPE=${2}
########################################################################
	cd ${Current_DATA_DIR}
	local DIRLIST="$( find -maxdepth 1 -name "*.${FILE_TYPE}" | sort -n | xargs)"
	
	local k=0
	for FILE_DIR in ${DIRLIST}
	do
		local _DIRLIST[k]="${Current_DATA_DIR}/${FILE_DIR: 2}"
		echo "Saving $k -th DIR as ${_DIRLIST[k]}"
		k=`expr $k + 1`
	done
	
	local Venn_Input_NUM=${#_DIRLIST[*]}
	echo "In total, ${Venn_Input_NUM} files are input."
	
	
	echo "Merge All files into one file!"
	echo "cat *.${FILE_TYPE} | cut -f 1,2,3 | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o count > union_all_0.txt"
	cat *.${FILE_TYPE} | cut -f 1,2,3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > union_all_0.txt
	
	
	if [ ${Venn_Input_NUM} == 2 ]
	then
		local k=0
		for FILE_DIR in ${DIRLIST}
		do
			local __DIR[k]="${FILE_DIR: 2}"
			local j=$(expr $k + 1)
			echo "bedtools intersect -c -a union_all_${k}.txt -b ${__DIR[k]} > union_all_${j}.txt"
			bedtools intersect -c -a union_all_${k}.txt -b ${__DIR[k]} > union_all_${j}.txt
			rm union_all_${k}.txt
			local k=$(expr $k + 1)
		done
		python ${Tools_DIR}/Venn_Diagram_plot2.py "union_all_${j}.txt" "${__DIR[0]::-4}" "${__DIR[1]::-4}"
	fi
	
	
	if [ ${Venn_Input_NUM} == 3 ]
	then
		local k=0
		for FILE_DIR in ${DIRLIST}
		do
			local __DIR[k]="${FILE_DIR: 2}"
			local j=$(expr $k + 1)
			echo "bedtools intersect -c -a union_all_${k}.txt -b ${__DIR[k]} > union_all_${j}.txt"
			bedtools intersect -c -a union_all_${k}.txt -b ${__DIR[k]} > union_all_${j}.txt
			rm union_all_${k}.txt
			local k=$(expr $k + 1)
		done
		python ${Tools_DIR}/Venn_Diagram_plot3.py "union_all_${j}.txt" "${__DIR[0]::-4}" "${__DIR[1]::-4}" "${__DIR[2]::-4}"
		rm union_all_${j}.txt
	fi
	#### Manually 
	#bedtools intersect -c -a union_all.${FILE_TYPE} -b TLE3-Tfh_peaks.bed | bedtools intersect -c -a stdin -b TLE3-Th1_peaks.bed | bedtools intersect -c -a stdin -b TLE3-CD4_peaks.bed > combine_sorted_count.bed 
	
	}

RUN_Bedtools_Intersect(){
	#### usage RUN_Bedtools_Intersect $1 $2
	CHECK_arguments $# 2
	local FILE_TYPE='bed'
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

RUN_Island_Filtered_Reads(){
	#### usage RUN_RPKM $1 $2
	CHECK_arguments $# 2
	local File_Name=${1}
	local FILE_TYPE=${2}
	
	local Output_Path="${__EXE_PATH}/island_filtered_reads"
	DIR_CHECK_CREATE ${Output_Path}
	
	echo "Entering the processing directory"
	cd ${__RAW_DATA_PATH_DIR}
	echo "From Reads bed file to get islandfiltered reads."
	echo "$(find -name "*${File_Name}*.${FILE_TYPE}" | sort -n | xargs)"
	local Input_A1="$( find -name "*${File_Name}*.${FILE_TYPE}" | sort -n | xargs)"
	local Input_A1="${__RAW_DATA_PATH_DIR}/${Input_A1: 2}"
	
	local Gene_list_folder=~/cloud_research/PengGroup/XLi/Annotation/MM9
	echo "Find islands from ${Gene_list_folder}"
	
	cd ${Gene_list_folder}
	echo "$( find -name "*.bed" | sort -n | xargs)"
	#local Input_Gene_Lists="$( find -name "Union_WT-na_dKO-na.bed" | sort -n | xargs)"
	local Input_Gene_Lists=( "./1070327_mm9_simple_repeat_and_Satellite.bed")
	
	
	for Input_Gene in ${Input_Gene_Lists}
	do
		#local OUTPUT_NAME="island_filtered_reads_${1}_${Input_Gene: 2: -4}"
		local OUTPUT_NAME="${1}_${Input_Gene: 2: -4}"
		Input_Gene="${Gene_list_folder}/${Input_Gene: 2}"
		echo "bedtools intersect -v -a ${Input_A1} -b ${Input_Gene} > ${Output_Path}/${OUTPUT_NAME}_repeat_removed.${FILE_TYPE}"
		bedtools intersect -v -a ${Input_A1} -b ${Input_Gene} > ${Output_Path}/${OUTPUT_NAME}_repeat_removed.${FILE_TYPE}
	done

}

RUN_RPKM(){
	#### usage RUN_RPKM $1 $2
	CHECK_arguments $# 2
	local File_Name=${1}
	local SPECIES=${2}
	local PATH_python_tools=~/cloud_research/PengGroup/XLi/Python_tools/read_RPKM.py
	
	
	
	echo "Entering the Raw Input directory"

	cd ${__RAW_DATA_PATH_DIR}

	echo "$(find -name "${File_Name}*.bedpe" | sort -n | xargs)"
	local Input_B1="$( find -name "*${File_Name}*.bedpe" | sort -n | xargs)"
	local Input_B1="${__RAW_DATA_PATH_DIR}/${Input_B1: 2}"
	local Input_Size=$(wc -l ${Input_B1} | cut -d' ' -f1)
	
	cd ${__EXE_PATH}
	echo "From Reads bed file to calculate reads count and its RPKM."
	case ${SPECIES} in
	"mm10")
	echo "Reference SPECIES is ${SPECIES}"
	local Gene_list_folder=~/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm10
	;;
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local Gene_list_folder=~/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm9
	local Gene_list_folder=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq/MACS2_Results/bampe_combined_repplicates
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
	
	#cd ${Gene_list_folder}
	#echo "$( find -name "*.${FILE_TYPE}" | sort -n | xargs)"
	#local Input_A1_Lists="$( find -name "*.${FILE_TYPE}" | sort -n | xargs)"
	
	local Input_A1_Lists=(
	Union_WT-na_dKO-na_1934_CTCF_top_2k_plus.bed
	)
	
	
	for Input_A1 in ${Input_A1_Lists[*]}
	do
		local OUTPUT_NAME="read_count_${1}_${Input_A1:: -4}"
		local Input_A1="${Gene_list_folder}/${Input_A1}"
		echo "cut -f 1,2,3,4 ${Input_A1} | bedtools intersect -c -a stdin -b ${Input_B1} > ${__EXE_PATH}/${OUTPUT_NAME}.bed"
		cut -f 1,2,3,4 ${Input_A1} | bedtools intersect -c -a stdin -b ${Input_B1} > ${__EXE_PATH}/${OUTPUT_NAME}.bed
		
		### A high memory efficient " -sorted " model requires both input files to be sorted. 
		#bedtools intersect -c -a ${Input_A1} -b ${Input_B1} > ${Output_Path}/${OUTPUT_NAME}.bed -sorted
		echo "python ${PATH_python_tools} ${OUTPUT_NAME} ${__EXE_PATH} ${Input_Size}"
		python ${PATH_python_tools} ${OUTPUT_NAME} ${__EXE_PATH} ${Input_Size}
	done

}

RUN_read_filtering_Bamtohdf5(){
####Usage: RUN_Bamtohdf5 $1(Input Parent folder path) $2(condition name) $3(Input subfolder) $4(Resolution) ...
	local SPECIES=mm9
	local MM9_GENOME_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM9/genome.fa
	local hg38_GENOME_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/hg38_UCSC/genome.fa
	local EXEDIR="~/Software/HiC_lib"
	local INPUT_Parameters=$@
	local INPUT_NUMS=$#
	local RESOLUTION=$4
	echo "Input ${INPUT_NUMS} parameters, which are: $*"
	
	local GENOME_PATH=${MM9_GENOME_PATH}
	
	for (( i = 3; i <= $(expr ${INPUT_NUMS} - 1); i++ ))
	do 
		local INPUT_R1=${1}/${2}/${3}/${2}_R1.bam
		local INPUT_R2=${1}/${2}/${3}/${2}_R2.bam
	
		local OUTPUT=${1}/${2}/${3}/f_hdf5
		DIR_CHECK_CREATE ${OUTPUT}
		
		
		echo "python run_parse_sam.py -g ${GENOME_PATH} -1 ${INPUT_R1} -2 ${INPUT_R2} -o ${OUTPUT}/${2}_mapped_reads.hdf5"
		#python ${EXEDIR}/run_parse_sam.py -g ${GENOME_PATH} -1 ${INPUT_R1} -2 ${INPUT_R2} -o ${OUTPUT}/${2}_mapped_reads.hdf5
		
		echo "python run_merge_pairs_and_rm_duplicates.py -g ${GENOME_PATH} -i ${OUTPUT} -o ${OUTPUT} -n ${2}"
		#python ${EXEDIR}/run_merge_pairs_and_rm_duplicates.py -g ${GENOME_PATH} -i ${OUTPUT} -n ${2} -o ${OUTPUT}
		
		echo "python run_ICE.py -g ${GENOME_PATH} -n ${2} -o ${OUTPUT} -r ${RESOLUTION}"
		#python ${EXEDIR}/run_ICE.py -g ${GENOME_PATH} -n ${2} -o ${OUTPUT} -r ${RESOLUTION}
		
		echo "		java -Xmx16g -jar ${EXEDIR}/juicer_tools_0.7.0.jar pre -n ${OUTPUT}/Juicebox/${2}_Juicebox_input.txt.gz ${OUTPUT}/Juicebox/${2}_Juicebox.hic $SPECIES"
		java -Xmx16g -jar ${EXEDIR}/juicer_tools_0.7.0.jar pre -n ${OUTPUT}/Juicebox/${2}_Juicebox_input.txt.gz ${OUTPUT}/Juicebox/${2}_Juicebox.hic $SPECIES
		
		echo "Finished Condition: $2"
		echo ""
		
		
		
	done
	
	
	
	
	
	
	
	echo "RUN_Bamtohdf5 is completed!"
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

RUN_SELECT_SEQ_SAM(){
	#### Usage: RUN_SELECT_SEQ_SAM $INPUT_SAM_FOLDER ${INPUT_NAME} ${NAME_EXTENSION} (Such as _unaligned for xxx_unaligned)
	#### Given a sam file DIR
	local INPUT_NAME=${1}
	local NAME_EXTENSION=${3}
	local INPUT_SAM_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
	local OUTPUT_RNAME=${2}
	CHECK_arguments $# 3
	
	echo "Entering one Library of SAM_DATA_FOLDER: $INPUT_SAM_FOLDER"
	cd ${INPUT_SAM_FOLDER}
	local INPUT_SAM=${1}${3}.sam
	local SORT_BAM=${1}${3}\_sorted.bam
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

RUN_REMOVE_SEQ(){
## USAGE:	RUN_REMOVE_SEQ ${__INPUT_SAMPLE_DIR_List[i]}
#### Check a specific seq in fastq.gz
####	1.Raw_data 2.Seq
local bbmap=${Tools_DIR}/bbmap
local REMOVE_SEQ=~/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/DNase_seq/Adapter_Remove/TruSeq_Adapter_Index_5.fq
local OUTPUT_NAME=${__RAW_DATA_PATH_DIR}/Adapter_Removed/${1}
DIR_CHECK_CREATE ${OUTPUT_NAME}
local KERS=${2}
########################################################################
echo ""
local DIRLIST_gz="$( find -name "*.fastq.gz" | sort -n | xargs)"

	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	echo "Pair End."
	echo 
	${bbmap}/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stats.txt
	echo ""
	#### concordantly pair output
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
	# using un concordantly pair-ends do bowtie2 again.
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}
	
	else
	echo "Single End."
	echo "bowtie2 -p $THREADS --no-unal --non-deterministic -x $BOWTIEINDEXS -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S $OUTPUT_BOWTIE2"
	${bbmap}/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stats.txt
	fi

	echo ""




echo "nohup $bbmap/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stas.txt"
nohup $bbmap/bbduk.sh in1=${__FASTQ_DIR_R1[0]} in2=${__FASTQ_DIR_R2[0]} out1=$OUT/unmatched_R1.fastq.gz out2=$OUT/unmatched_R2.fastq.gz outm1=$OUT/matched_R1.fastq.gz outm2=$OUT/matched2_R2.fastq.gz ref=$REMOVE_SEQ k=$KERS hdist=0 stats=stats.txt

########################################################################
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
	local Q_NAME=HWI
	
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

#RUN_Sam2Wig

RUN_bed2fastq(){
#### Usage: RUN_bed2fastq $1 $2($1 = Input_Name, $2 = FRAGMENT_SIZE) 
CHECK_arguments $# 2

local INPUT_PATH=${__RAW_DATA_PATH_DIR}
local OUT_PATH=${__EXE_PATH}/Associated_fastq
local INPUT_NAME=${1}
local SPECIES=${2}
########################################################################
DIR_CHECK_CREATE ${OUT_PATH}

### Find Input File
local IN_FOLDER=${__RAW_DATA_PATH_DIR}
cd ${IN_FOLDER}
local IN_FILES="$( find -name "*${INPUT_NAME}*.bed" | sort -n | xargs)"

local k=0
for in_files in ${IN_FILES}
do
	IN_FILES[k]="${IN_FOLDER}/${in_files: 2}"
	k=`expr $k + 1`
done
echo "Saving Bed INPUT File as ${IN_FILES[*]}"
### Find Input File

case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local FA_SEQUENCE=~/cloud_research/PengGroup/XLi/Annotation/MM9/mm9.fa
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local FA_SEQUENCE=~/cloud_research/PengGroup/XLi/Annotation/MM10/mm10.fa
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local FA_SEQUENCE=~/cloud_research/PengGroup/XLi/Annotation/HG19/hg19.fa
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local FA_SEQUENCE=~/cloud_research/PengGroup/XLi/Annotation/HG38/hg38.fa
	;;
	"dm6") 
	echo "Reference SPECIES is ${SPECIES} Not Setup Yet!"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac

### extension
local yesno='yes'
#awk '$2-=100' OFS='\t' ${input_file} | awk '$3+=100' OFS='\t' > ${input_file}_ex
#awk -F',' 'BEGIN {OFS='\t'} {$4 = "\tUnion_peaks_"(NR); print}' 28459_Union_peaks.bed >28459_Union_peaks_x.bed
cd ${OUT_PATH}
for input_file in ${IN_FILES[*]}
do
	echo "bedtools getfasta -name -fi ${FA_SEQUENCE} -bed ${input_file} -fo ${input_file::-4}.fastq"
	bedtools getfasta -fi ${FA_SEQUENCE} -bed ${input_file} -fo ${input_file::-4}.fastq -name
done



	}

RUN_bam2bedpe(){
	#### Usage: RUN_Sam2Wig $1 $2 ($1 = Input_Name, $2 = FRAGMENT_SIZE)
CHECK_arguments $# 2
local INPUT_PATH=${__RAW_DATA_PATH_DIR}/${1}/bowtie2_results/${1}
local OUT_PATH=${__EXE_PATH}/${1}
local FRAGMENT_SIZE=${2}
########################################################################

local SICER="/home/lxiang/Software/SICER1.1/SICER"
local EXEDIR="${SICER}/extra/tools/wiggle"
local WINDOW_SIZE=200
local SPECIES=mm10


#samtools view ${INPUT_PATH}.sam -Sb | bamToBed -i stdin > ${INPUT_PATH}.bed

bamToBed -i ${INPUT_PATH}.bam > ${INPUT_PATH}.bed

sh $EXEDIR/bed2wig.sh ${__RAW_DATA_PATH_DIR}/${1}/bowtie2_results ${1} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	
	}

RUN_Sam2Wig(){
	#### Usage: RUN_Sam2Wig $1 $2 ($1 = Input_Name, $2 = FRAGMENT_SIZE)
CHECK_arguments $# 2
local INPUT_PATH=${__RAW_DATA_PATH_DIR}/${1}/bowtie2_results/${1}
local OUT_PATH=${__EXE_PATH}/${1}
local FRAGMENT_SIZE=${2}
########################################################################

local SICER="/home/lxiang/Software/SICER1.1/SICER"
local EXEDIR="${SICER}/extra/tools/wiggle"
local WINDOW_SIZE=200
local SPECIES=mm10


#samtools view ${INPUT_PATH}.sam -Sb | bamToBed -i stdin > ${INPUT_PATH}.bed

bamToBed -i ${INPUT_PATH}.bam > ${INPUT_PATH}.bed

sh $EXEDIR/bed2wig.sh ${__RAW_DATA_PATH_DIR}/${1}/bowtie2_results ${1} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	
	}

RUN_BED2WIG(){
	#### Given a folder, change their all bed file to wig.
	#### RUN_BED2WIG $1 $2
	echo "RUN_BED2WIG"
	CHECK_arguments $# 2

	### Operation PARAMETERS Setting
	local INPUT_NAME=${1}
	local SPECIES=${2}
	
	local INPUT_DIR=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}
	
	#local SPECIES=hg19
	local WINDOW_SIZE=200
	local FRAGMENT_SIZE=50
	########################################################################
	########################################################################
	local EXEDIR=/opt/tools
	local SICER_DIR=${EXEDIR}/SICER1.1/SICER
	local Wiggle=${SICER_DIR}/extra/tools/wiggle
########################################################################
	local FILE_NAME=${__EXE_PATH}/${INPUT_NAME}.wig
	if [ ! -f ${FILE_NAME} ];then
	echo "sh $Wiggle/bed2wig.sh ${INPUT_DIR} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}"
	sh ${Wiggle}/bed2wig.sh ${INPUT_DIR} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	fi

#### Clear bed file.
	if [ -f ${FILE_NAME} ];then
	echo "rm ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bed"
	#rm ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bed
	fi
	
	}

RUN_Bed2BigBed(){
	### A copy of RUN_Wig2BigWig
	###RUN_Bed2BigBed $1 $2 $3 $4 $5 
	#### ($1=INPUT_FOLDER $2=INPUT_NAME $3=Track_hub_label, $4 Species $5 Data Provider)
	#### convert .bed to .bigbed and create the track hubs profiles.
	echo "RUN_Bed2BigBed"
	CHECK_arguments $# 5
	
	local Bed_DIR=${1}
	local Tracks_NAME=${2}
	local Data_label=${3}
	local SPECIES=${4}
	local Data_provider=${5}
	local Tracks_NAME_Lable=${Tracks_NAME: 7:12} ##Skip out of Sample_ (7) and forward 12

	####INPUT
case ${SPECIES} in 
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm9.chrom.sizes"
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm10.chrom.sizes"
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/hg19.chrom.sizes"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
	
	cd ${Bed_DIR}
#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigBeds


	echo "${UCSC_DIR}/bedToBigBed ${Tracks_NAME}.bed ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigBeds/${Tracks_NAME}.bb"
	${UCSC_DIR}/bedToBigBed ${Tracks_NAME}.bed ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigBeds/${Tracks_NAME}.bb

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

	local filename=genomes.txt
	if [ ! -f $filename ];then
	echo "genome $SPECIES" >>$filename
	echo "trackDb $SPECIES/trackDb.txt" >>$filename
	fi
########################################################################

	DIR_CHECK_CREATE $OUTPUTDIR_tracks_hub/${SPECIES}	
	cd $OUTPUTDIR_tracks_hub/${SPECIES}
	
	local filename=trackDb.txt
	local bw_Url=https://s3.amazonaws.com/xianglilab/tracks_hub/${Data_provider}/${Data_label}/BigBeds/${Tracks_NAME}.bb
	echo "track ${Tracks_NAME}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "type bigBed" >>$filename
	echo "visibility dense" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo "alwaysZero on" >>$filename
	echo "maxHeightPixels 100:16:8" >>$filename
	echo -e "bigDataUrl ${bw_Url} \n" >>$filename
########################################################################
########################################################################
#     Uncomleted! Apr.21.2018
########################################################################
########################################################################	
	}

RUN_Wig2BigWig (){
	###RUN_Wig2BigWig $1 $2 $3 $4 $5 
	#### ($1=INPUT_FOLDER $2=INPUT_NAME $3=Track_hub_label, $4 Species $5 Data Provider)
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $1 $2 $3
	CHECK_arguments $# 5
	
	local Wig_DIR=${1}
	local Tracks_NAME=${2}
	local Data_label=${3}
	local Data_provider=${5}
	local SPECIES=${4}
	local Tracks_NAME_Lable=${Tracks_NAME: 7:12} ##Skip out of Sample_ (7) and forward 12

	####INPUT
case ${SPECIES} in 
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm9.chrom.sizes"
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm10.chrom.sizes"
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/hg19.chrom.sizes"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac

	
	cd ${Wig_DIR}
	local Wig_File=$(find -name ${Tracks_NAME}*.wig | xargs)
#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigWigs


	echo "${UCSC_DIR}/wigToBigWig ${Wig_File: 2} ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	${UCSC_DIR}/wigToBigWig ${Wig_File: 2} ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw

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

	local filename=genomes.txt
	if [ ! -f $filename ];then
	echo "genome $SPECIES" >>$filename
	echo "trackDb $SPECIES/trackDb.txt" >>$filename
	fi
########################################################################

	DIR_CHECK_CREATE $OUTPUTDIR_tracks_hub/${SPECIES}	
	cd $OUTPUTDIR_tracks_hub/${SPECIES}
	
	local filename=trackDb.txt
	local bw_Url=https://s3.amazonaws.com/xianglilab/tracks_hub/${Data_provider}/${Data_label}/BigWigs/${Tracks_NAME}.bw
	echo "track ${Tracks_NAME}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "visibility full" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo "alwaysZero on" >>$filename
	echo "maxHeightPixels 100:24:8" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
	}

RUN_BedGraph2BigWig (){
	##RUN_BigGraph2BigWig   $INPUT_FOLDER $INPUT_NAME $Track_hub_label
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $Sample_Wig_NAME
	CHECK_arguments $# 5
	
	local Wig_DIR=${1}
	local Tracks_NAME=${2}
	local Data_label=${3}
	local Data_provider=${5}
	local SPECIES=${4}
	local Tracks_NAME_Lable=${Tracks_NAME: 7:12} ##Skip out of Sample_ (7) and forward 12
	####INPUT
case ${SPECIES} in 
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm9.chrom.sizes"
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/mm10.chrom.sizes"
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/hg19.chrom.sizes"
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local chrom_sizes="${UCSC_DIR}/genome_sizes/hg38.chrom.sizes"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac

	cd ${Wig_DIR}

#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigWigs
	
########################################################################
	
	if [ ! -f ${Tracks_NAME}_sorted.bdg ];then
	echo "sort -k1,1 -k2,2n ${Tracks_NAME}.bdg > ${Tracks_NAME}_sorted.bdg"
	sort -k1,1 -k2,2n ${Tracks_NAME}.bdg > ${Tracks_NAME}_sorted.bdg
	fi
	
	if [ -f ${Tracks_NAME}.bdg ];then
	rm ${Tracks_NAME}.bdg
	fi
	
	echo "$UCSC_DIR/bedGraphToBigWig ${Tracks_NAME}_sorted.bdg ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	$UCSC_DIR/bedGraphToBigWig ${Tracks_NAME}_sorted.bdg ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw

###	Creating hub.txt
	cd $OUTPUTDIR_tracks_hub
	local filename=hub.txt
	if [ ! -f $filename ];then
	echo "hub ${Data_label}" >>$filename
	echo "shortLabel ${Data_label}" >>$filename
	echo "longLabel ${Data_label}/${Data_provider}" >>$filename
	echo "genomesFile genomes.txt" >>$filename
	echo "email lux@gwu.edu" >>$filename
	echo "descriptionUrl descriptionUrl" >>$filename
	fi
########################################################################

###	Creating genomes.txt	
	local filename=genomes.txt
	if [ ! -f $filename ];then
	echo "genome $SPECIES" >>$filename
	echo "trackDb $SPECIES/trackDb.txt" >>$filename
	fi
########################################################################

	DIR_CHECK_CREATE $OUTPUTDIR_tracks_hub/${SPECIES}	
	cd $OUTPUTDIR_tracks_hub/${SPECIES}
	
	local filename=trackDb.txt
	local bw_Url=https://s3.amazonaws.com/xianglilab/tracks_hub/${Data_provider}/${Data_label}/BigWigs/${Tracks_NAME}.bw
	echo "track ${Tracks_NAME}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "visibility full" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo "alwaysZero on" >>$filename
	echo "maxHeightPixels 100:24:8" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
}

RUN_Reads_Profile(){
	### Usage: RUN_Reads_Profile $1 $2
	####
	
	
	CHECK_arguments $# 3
	#### TSS or TES or GeneBody
	local REGIONTYPE=${1}
	#### INPUT FILE NAME
	local INPUT_NAME=${2}
	local Gene_Type=${3} ### exon for retrosposon markers.  "transcript_region for IR gtf annotation file." 
	
	
	
	local FRAGMENTSIZE=100
	local UP_EXTENSION=5000
	local DOWN_EXTENSION=5000
	local WINDOWSIZE=100
	local RESOLUTION=10
	local NORMALIZATION=1.0
	local Genic_Partition=200
	
	
	#### Genelist for Profile
	local GENE_LIST_FOLDER=~/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/DNase_seq/MACS2_Results/bampe_combined_repplicates
	GENELIST_LIST=(
	Union_WT-na_dKO-na.bed
	)
	#### Genelist for Profile
	
	local EXE_PATH=~/cloud_research/PengGroup/XLi/Python_tools/generate_profile_TSS_GeneBody_TES_beta.py
	
	case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/2015_GTF/mm9_2015.gtf
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/2015_GTF/mm10_2015.gtf
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/2015_GTF/hg38_2015.gtf
	;;
	"Repeat_Masker")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Data/Paul/34bc/retrotransposons_marker/Choi-Data-S1.gtf
	;;
	
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
	
	
	echo "Entering the processing directory"
	cd ${__RAW_DATA_PATH_DIR}
	echo "From islandfilted Reads bed file to calculate Profile."
	echo "$(find -name "${INPUT_NAME}*-islandfiltered.bed" | sort -n | xargs)"
	local INPUT_FILE="$( find -name "${INPUT_NAME}*-islandfiltered.bed" | sort -n | xargs)"
	local INPUT_FILE="${__RAW_DATA_PATH_DIR}/${INPUT_FILE: 2}"
	
	local OUTPUTDIR=${__EXE_PATH}/${REGIONTYPE}_Profiles_Islandfiltered_Reads/${INPUT_NAME}
	DIR_CHECK_CREATE ${OUTPUTDIR}
	cd ${OUTPUTDIR}
	
	for GENELISTFILE in ${GENELIST_LIST[*]}
	do
	#local INPUT_FILE=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}/bowtie2_results/${INPUT_NAME}.bed
	#local INPUT_FILE=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bed
	echo "python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} -e ${Gene_Type} \
	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}"
	python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} -e ${Gene_Type} \
	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
	
	echo ""
	echo ""
	echo "done"
	#break
	done



	}

RUN_HomerTools(){
	### RUN_HomerTools $1 $2 $3 $4
	CHECK_arguments $# 2
	echo ""
	echo "RUN_HomerTools"
	local 
	local TRIM_TYPE=${1}
	local Restriction_Enzyme=${2}
	
	
	echo "3 <#> (trim this many bp off the 3' end of the sequence)"
	echo "-5 <#> (trim this many bp off the 5' end of the sequence)"
	for fastq_file in ${__FASTQ_DIR_R1[*]}
	do
		case ${TRIM_TYPE} in
		"dist")
		homerTools trim -3 0 -5 0 ${fastq_file}
		gzip ${fastq_file}.trimmed
		;;
		"Restriction_Enzyme")
		homerTools trim -3 ${Restriction_Enzyme} -mis 0 -matchStart 20 -min 20 ${fastq_file}
		gzip ${fastq_file}.trimmed
		;;
		"barcode")
		#homerTools barcodes
		;;
		"freq")
		#homerTools freq
		;;
		*)
		echo "ERR: Did Not Find Any Matched Reference...... Exit"
		exit
		;;
		esac
	done
	
	echo ""
	echo "One RUN_HomerTools is Completed!"
	}

RUN_Motif_Homer(){
	### RUN_Motif_Homer $1 $2 $3 $4
	CHECK_arguments $# 2
	echo ""
	echo "RUN_HomerTools"
	local 
	local TRIM_TYPE=${1}
	local Restriction_Enzyme=${2}
	
	
	echo "3 <#> (trim this many bp off the 3' end of the sequence)"
	echo "-5 <#> (trim this many bp off the 5' end of the sequence)"
	for fastq_file in ${__FASTQ_DIR_R1[*]}
	do
		case ${TRIM_TYPE} in
		"dist")
		homerTools trim -3 0 -5 0 ${fastq_file}
		gzip ${fastq_file}.trimmed
		;;
		"Restriction_Enzyme")
		homerTools trim -3 ${Restriction_Enzyme} -mis 0 -matchStart 20 -min 20 ${fastq_file}
		gzip ${fastq_file}.trimmed
		;;
		"barcode")
		#homerTools barcodes
		;;
		"freq")
		#homerTools freq
		;;
		*)
		echo "ERR: Did Not Find Any Matched Reference...... Exit"
		exit
		;;
		esac
	done
	
	echo ""
	echo "One RUN_HomerTools is Completed!"
	}
##END OF FUNDEMENTAL FUNCTIONS

######################################

########################################################################

######################################

##	MODOLES
########################################################################
########################################################################
## Alignor
RUN_BOWTIE2(){
	### RUN_BOWTIE2 $1 $2 $3 $4
	local INPUT_NAME=${1}
	local SPECIES=${2}
	local PROJECT_NAME=${3}
	local Data_Provider=${4}
	local just_align_yesno=${5}  ## This option is limit bowtie2 with only basic functions.
	
	CHECK_arguments $# 4
	echo ""
	echo "RUN_BOWTIE2"
	#### OUTPUT FORMAT
	local OUTPUT_BOWTIE2_FOLDER="${__EXE_PATH}/Bowtie2_Results/${INPUT_NAME}"
	DIR_CHECK_CREATE ${OUTPUT_BOWTIE2_FOLDER}
	
case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
	;;
	"dm6") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/Drosophila_Melanogaster/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
	;;
	"BOWTIEINDEXS_mm10_pMXs_combo")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Raw_Data/Paul/34bc/Bowtie2_Indexes/MM10_pMXs_combo/mm10_pMXs_combo
	;;
	
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac

########################################################################
	
	echo "#5’ (left)  3’ (right) end of each read "
    
	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
		echo "Pair End Mode"
		local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
		echo "bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") --trim3 0 --trim5 0 -S ${OUTPUT_BOWTIE2}" #--un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz 
		bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") --trim3 0 --trim5 0 -S ${OUTPUT_BOWTIE2} #--un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz 
		echo ""
		#### concordantly pair output
		
		#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
		#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
		# using un concordantly pair-ends do bowtie2 again.
		#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}"
		#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}
		
		###sam2bam   ##For MACS2
		if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
		echo "samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam"
		samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam
		fi
	else
		echo "Single End Mode."
		### cite ${Pair_Type} from Function PRE_READS_DIR
		local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
		echo "bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2} --trim3 0 --trim5 0"
		bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2} --trim3 0 --trim5 0
		
		### SINGLE END SAM TO BAM
		if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
		echo "samtools view -b ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam"
		samtools view -b ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam 
		fi
	fi
	
########################################################################
### Then clear sam file.
	if [ -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
		echo "rm ${OUTPUT_BOWTIE2}"
		rm ${OUTPUT_BOWTIE2}
	fi
########################################################################
	local just_align_yesno=${5}
	case ${just_align_yesno} in
	"yes")
	echo "Make Align Stop here, just need alignment!"
	echo ""
	return 0 ;;
	"no")
	echo "Continue!" ;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit;;
	esac
########################################################################	
	
	
## Remove Redundancy by Picard
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam ];then
		#echo "samtools sort -n -l 1 -o ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam "
		#samtools sort -n -l 1 -o ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam 

		echo "picard MarkDuplicates I=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam O=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam \
		M=${OUTPUT_BOWTIE2_FOLDER}/Marked_dup_metrics.txt REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true ASSUME_SORT_ORDER=queryname"
		picard MarkDuplicates I=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam O=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam \
		M=${OUTPUT_BOWTIE2_FOLDER}/Marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=queryname
		
		echo "rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam" 
		rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam
		
		#echo "rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam"
		#rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam
	fi

	


###bam2bed
	
	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	### In here, bedpe just represents the bed format of pair end reads.
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bedpe ];then
	echo "bamToBed -bedpe -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | cut -f 1,2,6,7 | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bedpe"
	bamToBed -bedpe -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | cut -f 1,2,6,7 | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bedpe
	
	echo "RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	fi
	else
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bed ];then
	echo "bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed"
	bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bed
	
	echo "RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	fi
	fi

	echo ""
	echo "One Bowtie2 is Completed!"
	}

RUN_TOPHAT(){
#### Usage: RUN_TOPHAT $1 $2 $3 $4

echo "RUN_TOPHAT_ANALYSIS"
CHECK_arguments $# 4

### Operation PARAMETERS Setting
local INPUT_NAME=${1}
local PROJECT_NAME=${2}
local SPECIES=${3}
local Data_Provider=${4}

#### OUTPUT FORMAT
local OUTPUT_TOPHAT_FOLDER="${__EXE_PATH}/Tophat_Results/${INPUT_NAME}"
DIR_CHECK_CREATE ${OUTPUT_TOPHAT_FOLDER}

local WINDOW_SIZE=200
local FRAGMENT_SIZE=50
########################################################################
case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/raw_mm9_genes.gtf
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
	;;
	"mm10")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/mm10_genes_raw.gtf
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/raw_hg19_genes.gtf
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
	;;
	"dm6") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/Drosophila_Melanogaster/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac


########################################################################


local EXEDIR=/opt/tools
local SICER_DIR=${EXEDIR}/SICER1.1/SICER
local Wiggle=${SICER_DIR}/extra/tools/wiggle

# Tang

#
########################################################################


#### Decide single end or pair ends mode
#### NOW it is only compatible with single file. Not with pieces files.
	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	echo "Pair End Mode"
	echo "tophat -p $THREADS -o ${OUTPUT_TOPHAT_FOLDER} --GTF ${GTFFILE} ${BOWTIEINDEXS} $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",")"
	tophat -p ${THREADS} -o ${OUTPUT_TOPHAT_FOLDER} --GTF ${GTFFILE} ${BOWTIEINDEXS} $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",")
	else
	echo "Single End Mode."
	echo "tophat -p $THREADS -o ${OUTPUT_TOPHAT_FOLDER} --GTF ${GTFFILE} ${BOWTIEINDEXS} $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",")"
	tophat -p ${THREADS} -o ${OUTPUT_TOPHAT_FOLDER} --GTF ${GTFFILE} ${BOWTIEINDEXS} $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",")
	fi

	echo "mv ${OUTPUT_TOPHAT_FOLDER}/accepted_hits.bam ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam"
	mv ${OUTPUT_TOPHAT_FOLDER}/accepted_hits.bam ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam
	
	echo "bedtools bamToBed -i ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam | sort -k1,1 -k2,2n > ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed"
	bedtools bamToBed -i ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam | sort -k1,1 -k2,2n > ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed
	
	echo "sh ${Wiggle}/bed2wig.sh ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}"
	sh ${Wiggle}/bed2wig.sh ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	
	### Then clear bed file.
	if [ -f ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed ];then
	echo "rm ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed"
	rm ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed
	fi

	
	echo "RUN_Wig2BigWig ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME}-W${WINDOW_SIZE}-RPKM ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Wig2BigWig ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME}-W${WINDOW_SIZE}-RPKM ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	
	echo "One TopHat is completed."
	echo ""
	}

RUN_CELLRANGER(){
#### Usage: RUN_TOPHAT $1 $2 $3

echo "RUN_CELLRANGER_ANALYSIS"
CHECK_arguments $# 2

### Operation PARAMETERS Setting
local INPUT_NAME=${1}
local SPECIES=${2}
#local Data_Provider=${4}



#### INPUT FORMAT
local INPUT_CELLRANGER=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}
#### OUTPUT FORMAT
local OUTPUT_TOPHAT_FOLDER="${__EXE_PATH}/CELLRANGER/${INPUT_NAME}"
DIR_CHECK_CREATE ${OUTPUT_TOPHAT_FOLDER}
cd ${OUTPUT_TOPHAT_FOLDER}
########################################################################
case ${SPECIES} in
	"mm10")
	echo "Reference SPECIES is ${SPECIES}"
	local TRANSCRIPTOME=~/cloud_research/PengGroup/XLi/Annotation/10X_GENOMICS/Cell_Ranger_Reference/refdata-cellranger-mm10-2.1.0
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/hg19_genes.gtf
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
	;;
	"dm6") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/Drosophila_Melanogaster/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
########################################################################

###
	echo "cellranger count --id=${PROJECT_NAME} \
	--transcriptome=${TRANSCRIPTOME} \
	--sample=${INPUT_NAME: 7} \
	--fastqs=${INPUT_CELLRANGER} \ "
	#--localcores=8 \
	#--localmem=32
	## --expect-cells=1000
	
	
	
	cellranger count --id=${INPUT_NAME} \
	--transcriptome=${TRANSCRIPTOME} \
	--sample=${INPUT_NAME: 7} \
	--fastqs=${INPUT_CELLRANGER} \
	#--localcores=16 \
	#--localmem=32
	#--expect-cells=1000
	
	echo "One CELLRANGER is completed."
	echo ""
	}

RUN_HiC_Pro(){
	### RUN_HiC_Pro $1 $2 $3 $4
	# HiC-Pro -i INPUT -o OUTPUT -c CONFIG [-s ANALYSIS_STEP] [-p] [-h] [-v]
	# [-s|--step ANALYSIS_STEP] : run only a subset of the HiC-Pro workflow; if not specified the complete workflow is run
    #  mapping: perform reads alignment
    #  proc_hic: perform Hi-C filtering
    #  quality_checks: run Hi-C quality control plots
    #  build_contact_maps: build raw inter/intrachromosomal contact maps
    #  ice_norm: run ICE normalization on contact maps
	local INPUT_NAME=${1}
	local SPECIES=${2}
	local PROJECT_NAME=${3}
	local Data_Provider=${4}

	CHECK_arguments $# 4
	echo ""
	echo "RUN_HiC_Pro"

	
	singularity exec hicpro_latest_ubuntu.img HiC-Pro -p
	
	
	
	#### OUTPUT FORMAT
	local OUTPUT_BOWTIE2_FOLDER="${__EXE_PATH}/Bowtie2_Results/${INPUT_NAME}"
	DIR_CHECK_CREATE ${OUTPUT_BOWTIE2_FOLDER}
	
case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/MM10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/HG38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
	;;
	"dm6") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Annotation/Drosophila_Melanogaster/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
	;;
	"BOWTIEINDEXS_mm10_pMXs_combo")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS=~/cloud_research/PengGroup/XLi/Raw_Data/Paul/34bc/Bowtie2_Indexes/MM10_pMXs_combo/mm10_pMXs_combo
	;;
	
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
########################################################################
	
	echo "#5’ (left)  3’ (right) end of each read "
    
	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
		echo "Pair End Mode"
		local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
		echo "bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") --trim3 1000 --trim5 0 -S ${OUTPUT_BOWTIE2}" #--un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz 
		#bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") --trim3 100 --trim5 0 -S ${OUTPUT_BOWTIE2} #--un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz 
		echo ""
		#### concordantly pair output
		
		#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
		#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
		# using un concordantly pair-ends do bowtie2 again.
		#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}"
		#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}
		
		###sam2bam   ##For MACS2
		if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
		echo "samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam"
		#samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam
		fi
	else
		echo "Single End Mode."
		### cite ${Pair_Type} from Function PRE_READS_DIR
		local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
		echo "bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2} --trim3 100 --trim5 5"
		bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2} --trim3 100 --trim5 5
		
		### SINGLE END SAM TO BAM
		if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
		echo "samtools view -b ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam"
		samtools view -b ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam 
		fi
	fi
## Remove Redundancy by Picard
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam ];then
	
	#echo "samtools sort -n -l 1 -o ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam "
	#samtools sort -n -l 1 -o ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam 
	
	echo "rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam" 
	#rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam
	
	echo "picard MarkDuplicates I=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam O=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam \
	M=${OUTPUT_BOWTIE2_FOLDER}/Marked_dup_metrics.txt REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true ASSUME_SORT_ORDER=queryname"
	picard MarkDuplicates I=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam O=${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam \
	M=${OUTPUT_BOWTIE2_FOLDER}/Marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=queryname
	
	#echo "rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam"
	#rm ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam
	fi

### Then clear sam file.
	if [ -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
	echo "rm ${OUTPUT_BOWTIE2}"
	rm ${OUTPUT_BOWTIE2}
	fi
	
########################################################################
########################################################################
	local yesno='no'
	case ${yesno} in
	"yes")
	echo "Make Align Stop here!"
	echo ""
	return 0 ;;
	"no")
	echo "Continue!" ;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit;;
	esac
########################################################################

###bam2bed
	
	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	### In here, bedpe just represents the bed format of pair end reads.
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bedpe ];then
	echo "bamToBed -bedpe -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | cut -f 1,2,6,7 | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bedpe"
	bamToBed -bedpe -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | cut -f 1,2,6,7 | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bedpe
	
	echo "RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	fi
	else
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed ];then
	echo "bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed"
	bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_Dup_Removed.bam | sort -k1,1 -k2,2n > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed
	
	echo "RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Bed2BigBed ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	fi
	fi

	echo ""
	echo "One Bowtie2 is Completed!"
	}


## Peaks Calling.
RUN_SICER(){
#### Usage: RUN_SICER $1 $2 $3 ($1 is input for SICER. $2 is the CONTRO Library, $3 is the Gap set.) 
echo "RUN_SICER"
CHECK_arguments $# 3
local EXEDIR="${Tools_DIR}/SICER1.1/SICER"
local FRAGMENT_SIZE=100
local REDUNDANCY=1
local WINDOW_SIZE=200
local INPUT_NAME=${1}
local INPUI_CON=${2}
local GAP_SET=${3}
local EFFECTIVEGENOME=0.85 
local FDR=0.05


local IN_SICER_FOLDER=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}
local IN_SICER_FILES=${INPUT_NAME}.bed

local CONTRO_SICER_DIR=${__RAW_DATA_PATH_DIR}/${INPUI_CON}
local CONTRO_SICER_FILE=${INPUI_CON}.bed

local OUT_SICER_FOLDER=${__EXE_PATH}/SICER_Results/${INPUT_NAME}
DIR_CHECK_CREATE ${OUT_SICER_FOLDER}

cd ${OUT_SICER_FOLDER}

echo "bash ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} mm9 ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log"
bash ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} "mm9" ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log



RUN_Wig2BigWig ${OUT_SICER_FOLDER} ${INPUT_NAME}-W200-normalized "Ezh2_ChIP_seq" "mm9" "Haihui"
RUN_Wig2BigWig ${OUT_SICER_FOLDER} ${INPUI_CON}-W200-normalized "Ezh2_ChIP_seq" "mm9" "Haihui"


	}

RUN_MACS2(){
#### Usage: RUN_MACS2 $1 $2 $3 $4 $5 ($1 is input for MACS. $2 is the CONTRO Library $3 is label) 
echo "RUN_MACS2"
CHECK_arguments $# 6


__RAW_DATA_PATH_DIR=${__EXE_PATH}

local INPUT_NAME=${1}
local INPUT_CON=${2}
local INPUT_LABEL=${3}
local SPECIES=${4}
local Data_Provider=${5}
local INPUT_TYPE=${6}
#local INPUT_LABEL=${INPUT_NAME: 7:3}
##Skip out of Sample_ (7), and forward with 4 more digits.

#local FDR=0.05
local p_value=0.00001

#### Default Input for MACS2 is bam format.
### Find Input File
local IN_FOLDER=${__RAW_DATA_PATH_DIR}
cd ${IN_FOLDER}
local IN_FILES="$( find -name "*${INPUT_NAME}*.bam" | sort -n | xargs)"
### Find Contro Files
local CONTRO_FOLDER=${__RAW_DATA_PATH_DIR}
cd ${CONTRO_FOLDER}
local CONTRO_FILE="$( find -name "*${INPUT_CON}*.bam" | sort -n | xargs)"

case ${INPUT_TYPE} in
	"bam")
	echo "bam Format"
	echo " "
	;;
	"bampe")
	echo "bampe Format"
	echo " "
	;;
	"bedpe")
	echo "Reference SPECIES is ${SPECIES}"
	local IN_FOLDER=${__RAW_DATA_PATH_DIR}
	cd ${IN_FOLDER}
	local IN_FILES="$( find -name "*${INPUT_NAME}*.bedpe" | sort -n | xargs)"
	### Find Contro Files
	local CONTRO_FOLDER=${__RAW_DATA_PATH_DIR}
	cd ${CONTRO_FOLDER}
	local CONTRO_FILE="$( find -name "*${INPUT_CON}*.bedpe" | sort -n | xargs)"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Type...... Exit"
	exit
	;;
esac

#### Search and save all input files PATH
local k=0
for in_files in ${IN_FILES}
do
	IN_FILES[k]="${IN_FOLDER}/${in_files: 2}"
	k=`expr $k + 1`
done
echo "Saving MACS2 INPUT File as ${IN_FILES[*]}"
### Find Input File

local k=0
for in_files in ${CONTRO_FILE}
do
	CONTRO_FILE[k]="${CONTRO_FOLDER}/${in_files: 2}"
	k=`expr $k + 1`
done
echo "Saving MACS2 Contro File as ${CONTRO_FILE[*]}"
### Find Contro Files

local OUT_FOLDER=${__EXE_PATH}/MACS2_Results/${INPUT_TYPE}/${INPUT_NAME}_vs_${INPUT_CON}
DIR_CHECK_CREATE ${OUT_FOLDER}
cd ${OUT_FOLDER}

case ${INPUT_TYPE} in
	"bam")
	echo "Reference SPECIES is ${SPECIES}"
	### --BAM
	echo "macs2 callpeak --format BAM -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME}_vs_${INPUT_CON} -B --SPMR -p ${p_value} " # --nomodel --shift -100 --extsize 200" 
	macs2 callpeak --format BAM -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME}_vs_${INPUT_CON} -B --SPMR -p ${p_value} # --nomodel --shift -100 --extsize 200
	;;
	"bampe")
	echo "Reference SPECIES is ${SPECIES}"
	### --BAMPE
	echo "macs2 callpeak --format BAMPE -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME}_vs_${INPUT_CON} -B -p ${p_value} --SPMR"  # --nomodel --shift -100 --extsize 200" 
	macs2 callpeak --format BAMPE -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME}_vs_${INPUT_CON} -B -p ${p_value} --SPMR # --nomodel --shift -100 --extsize 200
	;;
	"bedpe")
	echo "Reference SPECIES is ${SPECIES}"
	### --BEDPE  !!!!! BEDPE is a new formate of pair end reads, do not treat bed format as bedpe!!!!!
	echo "macs2 callpeak --format BEDPE -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME} -B --SPMR -p ${p_value}"
	macs2 callpeak --format BEDPE -t ${IN_FILES[*]} -c ${CONTRO_FILE[*]} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME}_vs_${INPUT_CON} -B --SPMR -p ${p_value}
	;;
	*)
	echo "ERR: Did Not Find Any Matched Type...... Exit"
	# for DNaseq To find enriched cutting sites such as some DNAse-Seq datasets. In this case, all 5' ends of sequenced reads should be extended in both direction to smooth the pileup signals. 
	# If the wanted smoothing window is 200bps, then use '--nomodel --shift -100 --extsize 200'.
	exit
	;;
esac


####      ${INPUT_NAME}_treat_pileup  for input

### Generating filtered Peaks from MACS2 output.
echo "python ${Python_Tools_DIR}/MACS2_peaks_filtering.py -i ${OUT_FOLDER} \
-n ${INPUT_NAME}_vs_${INPUT_CON}_peaks.xls -f 4.0 -p ${p_value} -q 0.05"

python ${Python_Tools_DIR}/MACS2_peaks_filtering.py -i ${OUT_FOLDER} \
-n ${INPUT_NAME}_vs_${INPUT_CON}_peaks.xls -f 4.0 -p ${p_value} -q 0.05

python ${Python_Tools_DIR}/MACS2_peaks_filtering.py -i ${OUT_FOLDER} \
-n ${INPUT_NAME}_vs_${INPUT_CON}_peaks.xls -f 0.0 -p ${p_value} -q 1.0

RUN_BedGraph2BigWig ${OUT_FOLDER} ${INPUT_NAME}_vs_${INPUT_CON}_treat_pileup ${INPUT_LABEL} ${SPECIES} ${Data_Provider}
RUN_BedGraph2BigWig ${OUT_FOLDER} ${INPUT_NAME}_vs_${INPUT_CON}_control_lambda ${INPUT_LABEL} ${SPECIES} ${Data_Provider}
}

RUN_MACS2_Diff(){  ###Need Updated
#### Usage: RUN_MACS2_Diff $1 $2 $3 $4 $5 $6 $threshold ($1 is input for con1. $2 is the CONTRO_1  $3 is Con2. $4 is the CONTRO_2 )
echo "RUN_MACS2 Diff"
CHECK_arguments $# 2
local EXEDIR="~/Software/python_tools/MACS2-2.1.1.20160309/bin"

local CON1_NAME=${1}
local CON2_NAME=${2}


########################################################################
### Find Input File
local IN_FOLDER=${__RAW_DATA_PATH_DIR}
cd ${IN_FOLDER}
local IN_FILES="$( find -name "*${CON1_NAME}*.bdg" | sort -n | xargs)"

local k=0
for in_files in ${IN_FILES}
do
	IN_FILES[k]="${IN_FOLDER}/${in_files: 2}"
	k=`expr $k + 1`
done
echo "Saving MACS2 Diff INPUT File as ${IN_FILES[*]}"
### Find Input File

### Find Contro Files
local CONTRO_FOLDER=${__RAW_DATA_PATH_DIR}
cd ${CONTRO_FOLDER}
local CONTRO_FILE="$( find -name "*${CON2_NAME}*.bdg" | sort -n | xargs)"

local k=0
for in_files in ${CONTRO_FILE}
do
	CONTRO_FILE[k]="${CONTRO_FOLDER}/${in_files: 2}"
	k=`expr $k + 1`
done
echo "Saving MACS2 Diff Contro File as ${CONTRO_FILE[*]}"
### Find Contro Files

########################################################################

########################################################################
local OUT_FOLDER=${__EXE_PATH}/MACS2_Diff_results/${CON1_NAME}_vs_${CON2_NAME}
DIR_CHECK_CREATE ${OUT_FOLDER}

#This technique allows for a variable to be assigned a value if another variable is either empty or is undefined

local Peaks_File="$( find -name "*${CON1_NAME}*_peaks.xls" | sort -n | xargs)"
local Size_Treat1=$(egrep "fragments after filtering in treatment" ${IN_FOLDER}/${Peaks_File: 2} | grep -o '[0-9]*')
local Size_Contro1=$(egrep "fragments after filtering in control" ${IN_FOLDER}/${Peaks_File: 2} | grep -o '[0-9]*')
local Size_Contro1="${Size_Contro1:-0}"


local Peaks_File="$( find -name "*${CON2_NAME}*_peaks.xls" | sort -n | xargs)"
local Size_Treat2=$(egrep "fragments after filtering in treatment" ${CONTRO_FOLDER}/${Peaks_File: 2} | grep -o '[0-9]*')
local Size_Contro2=$(egrep "fragments after filtering in control" ${CONTRO_FOLDER}/${Peaks_File: 2} | grep -o '[0-9]*')
local Size_Contro2="${Size_Contro2:-0}"

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

cd ${OUT_FOLDER}


echo "macs2 bdgdiff --t1 ${IN_FILES[1]} --c1 ${IN_FILES[0]} --t2 ${CONTRO_FILE[1]} \
--c2 ${CONTRO_FILE[0]} --d1 ${Size_Treat1} --d2 ${Size_Treat2} -g 60 -l 120 --o-prefix diff_${CON1_NAME}_vs_${CON2_NAME}"


macs2 bdgdiff --t1 ${IN_FILES[1]} --c1 ${IN_FILES[0]} --t2 ${CONTRO_FILE[1]} \
--c2 ${CONTRO_FILE[0]} --d1 ${Size_Treat1} --d2 ${Size_Treat2} -g 60 -l 120 --o-prefix diff_${CON1_NAME}_vs_${CON2_NAME}

}

RUN_CUFFDIFF(){
	### RUN_CUFFDIFF $1 $2
	####
CHECK_arguments $# 4
local INPUT_Args=("$@")
local SAMPLE_NUM=${#INPUT_Args[*]}
local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/mm10_genes.gtf
local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/mm9_genes.gtf
local OUTPUT_Cuffdiff=${__EXE_PATH}/Cuffdiff_Results
DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}

#### Preparing the .bam files for cuffdiff
for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do	
	cd ${__RAW_DATA_PATH_DIR}
	local DATA_BAM[$i]="$( find -name "${INPUT_Args[i]}.bam" | sort -n | xargs)"
	echo ""
	echo "CuffDiff Library"
	echo ${DATA_BAM[$i]}
	echo ""
done

### Four GROUP, each group has 3 DATA FILES, 
########################################################################
if [ ${SAMPLE_NUM} -eq "4" ]
then
	local A_Label="WT_Tfh_Day8"
	local B_Label="EKO_Tfh_Day8"
	DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}/${B_Label}_vs_${A_Label}
	echo "${INPUT_Args} Data files are loading..."
	echo "cuffdiff -q -o ${OUTPUT_Cuffdiff}/${B_Label}_vs_${A_Label} -p ${THREADS} -L "${A_Label}","${B_Label}" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]} ${DATA_BAM[2]},${DATA_BAM[3]}"
	cuffdiff -q -o ${OUTPUT_Cuffdiff}/${B_Label}_vs_${A_Label} -p ${THREADS} -L "${A_Label}","${B_Label}" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]} ${DATA_BAM[2]},${DATA_BAM[3]}
	echo ""
	echo ""
fi

########################################################################

echo "CuffDiff Completed!"
	}

RUN_Quant_IRI(){
	#### Usage: RUN_Quant_IRI 1.INPUT_FOLDER 2.${INPUT_SAMPLE_DIR_List[*]}
	INPUT_FOLDER=${1}
	INPUT_NAME=${2}
	echo "RUN_Quant_IRI"
	CHECK_arguments $# 2
	local MAPDIR=/home/lxiang/cloud_research/PengGroup/ZZeng/Annotation/mappability
	local MAPFILE=hg19_wgEncodeCrgMapabilityAlign50mer.bigWig
	
	echo "IRTools quant -q IRI -i ${INPUT_FOLDER}/${INPUT_NAME}.bam -p single -s fr-unstranded -e hg19 -u ${MAPDIR}/${MAPFILE} -n ${INPUT_NAME}"
	IRTools quant -q IRI -i ${INPUT_FOLDER}/${INPUT_NAME}.bam -p single -s fr-unstranded -e hg19 -u ${MAPDIR}/${MAPFILE} -n ${INPUT_NAME}
	echo "Quant for Library: ${INPUT_NAME} is finished."
	}

RUN_Peaks_Distribution_Analysis(){
	####RUN_Peaks_Distribution_Analysis $1
	CHECK_arguments $# 2
	
	local INPUT_NAME=${1}
	local INPUT_FILE_TYPE=${2}
	
	local EXEDIR=${Tools_DIR}/Python_tools

	local GTFFILE=~/cloud_research/PengGroup/XLi/Annotation/gtf_files/mm10_genes.gtf

### This parameter means that Promoter [-1k TSS, +1k TSS] Gene_body[+1k TSS. TES]
	local PROMOTER_UPSTREAM_EXTENSION=1001   # Upstream extension, the distance from TSS.
	local TSS_REGION_LENGTH=2000


	local INPUTFILE=${INPUT_NAME}.${INPUT_FILE_TYPE}
	local OUTPUTFILE=${INPUT_NAME}_distribution.txt

	echo "python ${EXEDIR}/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g ${GTFFILE} -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE"
	python ${EXEDIR}/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g ${GTFFILE} -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE
	echo ""
	echo ""

	}

RUN_HiC_Iterative_Mapping(){
#### Usage: RUN_HiC_Iterative_Mapping #1 #2
	CHECK_arguments $# 1
	echo "RUN HiC Iterative Mapping!"
########################################################################
	local EXEDIR="/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/version.beta"
	local OUTPUT_BAM_PATH="${__EXE_PATH}/$1/iterative_mapping_bam"
########################################################################
	local BOWTIE_PATH="/usr/local/bowtie2-2.2.6/bin/bowtie2"
	local BOWTIE_INDEX_PATH="/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
	local BOWTIE_INDEX_PATH="/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
	DIR_CHECK_CREATE ${EXEDIR} ${OUTPUT_BAM_PATH}
	
	local IM_PARAMETERS=(25 10 0 50) # min_seq_len, len_step, seq_start, seq_end
	
	# ${__FASTQ_DIR_R1[0]}
	
	
	
	for fastq_path in ${__FASTQ_DIR_R1}
	do
	echo "python $EXEDIR/run_iterative_mapping.py -b ${BOWTIE_PATH} -i ${BOWTIE_INDEX_PATH} -q ${fastq_path}\
	-o ${OUTPUT_BAM_PATH}/${__INPUT_SAMPLE_DIR_List[i]}_R1.bam -m ${IM_PARAMETERS[0]} -t ${IM_PARAMETERS[1]} -s ${IM_PARAMETERS[2]} -e ${IM_PARAMETERS[3]}"
	python $EXEDIR/run_iterative_mapping.py -b ${BOWTIE_PATH} -i ${BOWTIE_INDEX_PATH} -q ${fastq_path}\
	-o ${OUTPUT_BAM_PATH}/${__INPUT_SAMPLE_DIR_List[i]}_R1.bam -m ${IM_PARAMETERS[0]} -t ${IM_PARAMETERS[1]} -s ${IM_PARAMETERS[2]} -e ${IM_PARAMETERS[3]}
	done
	
	for fastq_path in ${__FASTQ_DIR_R2}
	do
	echo ""
	python $EXEDIR/run_iterative_mapping.py -b ${BOWTIE_PATH} -i ${BOWTIE_INDEX_PATH} -q ${fastq_path}\
	-o ${OUTPUT_BAM_PATH}/${__INPUT_SAMPLE_DIR_List[i]}_R2.bam -m ${IM_PARAMETERS[0]} -t ${IM_PARAMETERS[1]} -s ${IM_PARAMETERS[2]} -e ${IM_PARAMETERS[3]}
	done
	
	
	echo "Finished!..."
	
	}
## END OF MODULES

######################################

########################################################################

######################################

##BASIC FUNCTIONS
########################################################################
########################################################################

FUNC_TEST(){
	local AAA=$1
	local BBB=$2
	CHECK_arguments $# 2
	echo "function AAA is $AAA"
	read -p "Input a number from keyboard." Read_key
	#echo "Your input is: ${Read_key}"
	printf "Your input is: %s ${Read_key} \n"
} 

FUNC_BED_Sort(){
	CHECK_arguments $# 2
	Tracks_NAME=${1}
	File_Type=${2}
	
	if [ ! -f ${Tracks_NAME}_sorted.${File_Type} ];then
	echo "sort -k1,1 -k2,2n ${Tracks_NAME}.${File_Type} > ${Tracks_NAME}_sorted.${File_Type}"
	sort -k1,1 -k2,2n ${Tracks_NAME}.${File_Type} > ${Tracks_NAME}_sorted.${File_Type}
	
	echo "After sort, delete oringinal file"
	rm ${Tracks_NAME}.${File_Type}
	fi
	
	
	}

FUNC_CHOOSE_EMAIL_ALERT(){
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
		echo "No Email Alert Confirmed."
		return 1 
		break
		else
		echo "No Alert Confirmed."
		return 0
		break
		fi
	done
}

EMAIL_ME(){
	CHECK_arguments $# 2
	echo "Start at ${1} " | mail -s "Project: + ${2} Finished" lux@gwu.edu
	}

FUNC_Download (){
	CHECK_arguments $# 1
### Download Login information and the download directory.
#### -nH --cut-dirs=3   Skip 3 directory components.
	local Web_Address="http://dnacore454.healthcare.uiowa.edu/20180801-0173_Xue_Qiang10xlHvkfpSaYwqguDJdpweNNwvWQABoPGDdUkPXVGMZ/results/Project_Xue_Qiang10x/"
	local Directory_Skip_Num=4
	local USER="gec"
	local USER=""
	local PASSWORD="aardvark dryer rummage"
	local PASSWORD=""
	cd ${__RAW_DATA_PATH_DIR}
	local DOWNLOAD_TARGET_NAME=${1}/;
####If download file is a folder. IT MUST END WITH trailing slash  "/"
########################################################################



	local Down_dir=${__RAW_DATA_PATH_DIR}/${DOWNLOAD_TARGET_NAME}/;
	
	if [ -n "${USER}" -a -n "${PASSWORD}" ]
	then
	echo "wget --no-check-certificate -nv -r -c -nH --cut-dirs=${Directory_Skip_Num} --user=${USER} --password=${PASSWORD} --accept=gz --no-parent ${Web_Address}/${DOWNLOAD_TARGET_NAME}"
	wget --no-check-certificate -nv -r -c -nH --cut-dirs=${Directory_Skip_Num} --user="${USER}" --password="${PASSWORD}" --accept=gz --no-parent ${Web_Address}/${DOWNLOAD_TARGET_NAME}
	else
#### NO PASSCODE
	echo "wget --no-check-certificate -nv -r -c -nH --cut-dirs=${Directory_Skip_Num} --accept=gz --no-parent ${Web_Address}/${DOWNLOAD_TARGET_NAME}"
	wget --no-check-certificate -nv -r -c -nH --cut-dirs=${Directory_Skip_Num} --accept=gz --no-parent ${Web_Address}/${DOWNLOAD_TARGET_NAME}
	fi
	echo "${DOWNLOAD_TARGET_NAME} Downloaded Completed"
	echo ""

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

DIR_CHECK_CREATE(){
### Saving DIR Check and Create
#### Usage: DIR_CHECK_CREATE $@
	echo ""
	echo "DIR_CHECK_CREATE!"
	echo "Input Dir: $@"
	local Dirs=$@
	for solo_dir in $Dirs
	do
		if [ ! -d $solo_dir ];then
			echo "Dir check and create is $solo_dir"
			echo ""
			mkdir -p $solo_dir
		fi
	done
	}

FUNC_BACKUP_FOLDER (){
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

FUN_GZIP(){
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

FUNC_CLEAN(){
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

FUNC_CUT_Columns (){
	####RUN_CUT_Columns $Filename $start $end
	CHECK_arguments $# 3
	local File_type='bed'
	local File_Path=${__RAW_DATA_PATH_DIR}/${1}.${File_type}
	local start=${2}
	local end=${3}
	local delimiter="	"
	cut -f ${start}-${end} ${File_Path} > ${__RAW_DATA_PATH_DIR}/${1}_Columns_${start}_${end}.${File_type}
	echo ""
	echo ""
	}

FUNC_CUT_Rows (){
	####FUNC_CUT_Rows $Filename $start $end
	
	#### Normally sam file format, '1,24d' is removing all header.
	#CHECK_arguments $# 3
	local File_type='sam'
	local File_Path=${__RAW_DATA_PATH_DIR}/${1}.${File_type}
	local start=${2}
	local end=${3}
	sed '${start},${end}d' ${File_Path} > ${__RAW_DATA_PATH_DIR}/${1}_Row_${start}_${end}.${File_type}
	echo ""
	echo ""
}

FUNC_CUT_Rows_Redirect (){
	####FUNC_CUT_Rows $Filename $start $end
	
	#### Normally sam file format, '1,24d' is removing all header.
	#CHECK_arguments $# 3
	local File_type='sam'
	local File_Path=${__RAW_DATA_PATH_DIR}/${1}.${File_type}
	local start=${2}
	local end=${3}
	sed -e '${start},${end}d' ${File_Path} > ${__RAW_DATA_PATH_DIR}/${1}_Row_${start}_${end}.${File_type}
	echo "a"
	echo ""
}


FUNC_Max(){
## Usage: Larger_value = $(RUN_Max $num_a $num_b)
	CHECK_arguments $# 2
	if [ $1 -gt $2 ]
	then
		echo $1
	else
		echo $2
	fi
	}
##END OF BASIC FUNCTIONS
########################################################################
########################################################################

######################################

########################################################################

######################################

##### Following Line is very IMPORTANT 
#main "$@"
