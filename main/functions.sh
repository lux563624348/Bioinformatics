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

######################################

##	FUNDEMENTAL FUNCTIONS FOR ANALYSIS MODULES
PRE_READS_DIR(){
	### PRE_READS_DIR ${__INPUT_SAMPLE_DIR_List[i]} "fastq.gz"
	CHECK_arguments $# 2
	echo "Entering one Library of RAW_DATA_DIR: $__RAW_DATA_PATH_DIR/$1"
	echo "Searching file type as: $2"
	local Current_DATA_DIR="${__RAW_DATA_PATH_DIR}/$1"
	cd ${Current_DATA_DIR}
	
	
########################################################################
	#local DIRLIST_R1="$( find -name "*R1.fastq" | sort -n | xargs)"
	#local DIRLIST_R2="$( find -name "*R2.fastq" | sort -n | xargs)" #
	
	local DIRLIST_R1_gz="$( find -name "*R1*.$2" | sort -n | xargs)"
	local DIRLIST_R2_gz="$( find -name "*R2*.$2" | sort -n | xargs)"


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

RUN_COPY_AND_CHANGE_NAME(){
	### RUN_COPY_AND_CHANGE_NAME $1 $2 $3 ($1 is INPUT_NAME $2> OUTPUT_NAME $3 "fastq.gz")
	CHECK_arguments $# 3
	local INPUT_DATA_PATH="${__RAW_DATA_PATH_DIR}/${1}"
	local OUTPUT_DATA_PATH="${__EXE_PATH}/${2}"
	
	mkdir -p ${OUTPUT_DATA_PATH}
	
	cp ${INPUT_DATA_PATH}.${3} ${OUTPUT_DATA_PATH}/${2}.${3}
	echo "Finished Copying File: ${2}"
	
	
	
	}

RUN_FAST_QC (){
####	RUN_FAST_QC
####	Usage: RUN_FAST_QC $INPUT_DATA_DIR
########################################################################
### FASTQC Output DIR setting ...
	local Output_fastqc=${__EXE_PATH}/fastqc
	DIR_CHECK_CREATE ${Output_fastqc}
	cd ${__EXE_PATH}
	
	for fastq_file in ${__FASTQ_DIR_R1[*]}
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

RUN_BedToFasta(){
	#### usage RUN_BedToFasta $1 $2  ($1 is folder $2 is name)
	#CHECK_arguments $# 2
	
	####Fasta files:
	local MM9_FASTA="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM9/genome.fa"
	local MM10_FASTA="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM10/mm10.fa"
	####
	local Input_1=${1}/${2}
	
	echo "bedtools getfasta -fi ${MM9_FASTA} -bed ${Input_1}"
	bedtools getfasta -fo ${Input_1}.fa -fi ${MM9_FASTA} -bed ${Input_1}.bed
	}

RUN_Venn_Diagram(){
	#### Usage: RUN_Bedtools_Merge $1 $2
	:'
	For example:
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

	Because chr1	1	30 has overlap both in ‘A’ and ‘B’, then we says that this peak is common in A and B.
	'
	
	CHECK_arguments $# 2
	local Current_DATA_DIR=${1}
	local FILE_TYPE=${2}
########################################################################
	cd ${Current_DATA_DIR}
	local DIRLIST="$( find -name "*.${FILE_TYPE}" | sort -n | xargs)"
	
	echo "Merge All files into one file!"
	echo "cat *.${FILE_TYPE} | sort -k1,1 -k2,2n | bedtools merge -i stdin > union_all_0.txt"
	cat *.${FILE_TYPE} | sort -k1,1 -k2,2n | bedtools merge -i stdin > union_all_0.txt
	
	
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
	
	#### Manually 
	#bedtools intersect -c -a union_all.${FILE_TYPE} -b TLE3-Tfh_peaks.bed | bedtools intersect -c -a stdin -b TLE3-Th1_peaks.bed | bedtools intersect -c -a stdin -b TLE3-CD4_peaks.bed > combine_sorted_count.bed 
	
	python Venn_Diagram_plot.py "union_all_${j}.txt" "${__DIR[0]::-4}" "${__DIR[1]::-4}" "${__DIR[2]::-4}"
	
	
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

RUN_RPKM(){
	#### usage RUN_Bedtools_Intersect $1 $2
	#### usage RUN_Bedtools_Intersect $1 $2
	CHECK_arguments $# 1
	local FILE_TYPE='bed'
	local PATH_python_tools="/home/lxiang/cloud_research/PengGroup/XLi/Python_tools/read_RPKM.py"
	local Input_A1="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm9/gene_Genebody_unique.bed"
	local Input_B1=${__RAW_DATA_PATH_DIR}/${1}/Tophat_results/${1}.${FILE_TYPE}

	
	####################################################################
	echo ""
	
	
	Output_Path="${__EXE_PATH}/genes_RPKM"
	DIR_CHECK_CREATE ${Output_Path}
	
	
	OUTPUT_NAME="genes_read_count_${1}.bed"
	
	
	

	echo "bedtools intersect -c -a ${Input_A1} -b ${Input_B1} > ${Output_Path}/${OUTPUT_NAME}"
	bedtools intersect -c -a ${Input_A1} -b ${Input_B1} > ${Output_Path}/${OUTPUT_NAME}
	
	echo "python ${PATH_python_tools} ${OUTPUT_NAME} ${Output_Path}"
	python ${PATH_python_tools} ${OUTPUT_NAME} ${Output_Path}

}

RUN_read_filtering_Bamtohdf5(){
####Usage: RUN_Bamtohdf5 $1(Input Parent folder path) $2(condition name) $3(Input subfolder) $4(Resolution) ...
	local SPECIES=mm9
	local MM9_GENOME_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM9/genome.fa
	local hg38_GENOME_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/hg38_UCSC/genome.fa
	local EXEDIR="/home/lxiang/Software/HiC_lib"
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

	#local SPECIES=hg19
	local WINDOW_SIZE=200
	local FRAGMENT_SIZE=50
	local THREADS=4
	########################################################################
	#### mm9
	local GTFFILE=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2014-05-23-16-05-24/Genes/genes.gtf
	local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

	#### hg19
	#local GTFFILE=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf
	#local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

	########################################################################
	local SICER_DIR=/home/lxiang/Software/SICER1.1/SICER
	local EXEDIR=$SICER_DIR/extra/tools/wiggle
########################################################################
	local FILE_NAME=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.wig
	if [ ! -f ${FILE_NAME} ];then
	echo "sh $EXEDIR/bed2wig.sh ${__EXE_PATH} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}"
	sh $EXEDIR/bed2wig.sh ${__EXE_PATH} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	fi
	
	}

RUN_Wig2BigWig (){
	###RUN_Wig2BigWig $1 $2 $3 $4 $5 
	#### ($1=INPUT_FOLDER $2=INPUT_NAME $3=Track_hub_label, $4 Species $5 Data Provider)
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $1 $2 $3
	CHECK_arguments $# 5
	
	local Data_provider=${5}
	local SPECIES=${4}
	
	local Wig_DIR=${1}
	local Tracks_NAME=${2}
	local Data_label=${3}
	local Tracks_NAME_Lable=${Tracks_NAME: 7:12} ##Skip out of Sample_ (7) and forward 12
	local UCSC_DIR=/home/lxiang/Software/UCSC
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
#### OUTPUT
	local OUTPUTDIR_tracks_hub=${__EXE_PATH}/tracks_hub/${Data_provider}/${Data_label}
	DIR_CHECK_CREATE ${OUTPUTDIR_tracks_hub}/BigWigs


	echo "$UCSC_DIR/wigToBigWig ${Tracks_NAME}.wig ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	$UCSC_DIR/wigToBigWig ${Tracks_NAME}.wig ${chrom_sizes} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw

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
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "visibility full" >>$filename
	echo "maxHeightPixels 128" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
	}

RUN_Bed2BigBed(){
	#### Given a folder, change their all bed file to BigBed.
	#### RUN_BED2BigBed $1 $2
	echo "RUN_BED2WIG"
	CHECK_arguments $# 2

	### Operation PARAMETERS Setting
	local INPUT_NAME=${1}
	local SPECIES=${2}

	local WINDOW_SIZE=200
	local FRAGMENT_SIZE=50
	local THREADS=4
	########################################################################
	#### mm9
	local GTFFILE=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2014-05-23-16-05-24/Genes/genes.gtf
	local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

	#### hg19
	#local GTFFILE=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf
	#local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

########################################################################
	local UCSC_DIR=/home/lxiang/Software/UCSC
	#local chrom_sizes_mm10=${UCSC_DIR}/genome_sizes/mm10.chrom.sizes
	local chrom_sizes_mm9=${UCSC_DIR}/genome_sizes/mm9.chrom.sizes
	#local chrom_sizes_hg19=${UCSC_DIR}/genome_sizes/hg19.chrom.sizes

########################################################################
	local FILE_NAME=${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bigBed
	
	if [ ! -f ${FILE_NAME} ];then
		echo "sort -kl,1 -k2,2n ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bed > ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}_sorted.bed"
		sort -kl,1 -k2,2n ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}.bed > ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}_sorted.bed
		
		
	
	fi
	
	echo "$UCSC_DIR/bedGraphToBigWig ${Tracks_NAME}_sorted.bdg $Reference_genome_sizes ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw"
	$UCSC_DIR/bedToBigBed ${__RAW_DATA_PATH_DIR}/${INPUT_NAME}_sorted.bed ${chrom_sizes_mm9} ${OUTPUTDIR_tracks_hub}/BigWigs/${Tracks_NAME}.bw
########################################################################
########################################################################
#     Uncomleted! Apr.21.2018
########################################################################
########################################################################	
	}

RUN_BigGraph2BigWig (){
	##RUN_BigGraph2BigWig   $INPUT_FOLDER $INPUT_NAME $Track_hub_label
#### convert .wig to .bigwig and create the track hubs profiles.
#### Usage: RUN_Wig2BigWig $Sample_Wig_NAME
	CHECK_arguments $# 5
	
	local Data_provider=${5}
	local SPECIES=${4}
	
	local Wig_DIR=${1}
	local Tracks_NAME=${2}
	local Data_label=${3}
	local Tracks_NAME_Lable=${Tracks_NAME: 7:12} ##Skip out of Sample_ (7) and forward 12
	local UCSC_DIR=/home/lxiang/Software/UCSC
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

	cd ${bdg_DIR}

	#local chrom_sizes_mm10=${UCSC_DIR}/genome_sizes/mm10.chrom.sizes
#	local chrom_sizes_mm9=${UCSC_DIR}/genome_sizes/mm9.chrom.sizes
	#local chrom_sizes_hg19=${UCSC_DIR}/genome_sizes/hg19.chrom.sizes


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
	echo "type bigWig" >>$filename
	echo "bigDataUrl ${bw_Url}" >>$filename
	echo "shortLabel ${Tracks_NAME: 0:18}" >>$filename
	echo "longLabel ${Tracks_NAME}" >>$filename
	echo "visibility full" >>$filename
	echo "maxHeightPixels 128" >>$filename
	echo "color 256,0,0" >>$filename
	echo "autoScale on" >>$filename
	echo -e "windowingFunction mean+whiskers \n" >>$filename
}

RUN_Reads_Profile(){
	### Usage: RUN_Reads_Profile $1 
	####
	
	
	CHECK_arguments $# 1
	
	local REGIONTYPE=${1}
	
	local FRAGMENTSIZE=0
	local UP_EXTENSION=5000
	local DOWN_EXTENSION=5000
	local WINDOWSIZE=100
	local RESOLUTION=1
	local FOLDER_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Hdac/Treg_HistoneMarks_ChIPSeq_Sep2016/SICER/WT_Treg_D20_H3K9Ac
	local GENE_LIST_FOLDER=${FOLDER_PATH}
	
	
	local EXE_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Python_tools/generate_profile_around_locations.py
	local GTFFILE=/home/lxiang/cloud_research/PengGroup/ZZeng/Annotation/gtf_files/mm9_genes.gtf
	
	
	local OUTPUTDIR=${__EXE_PATH}/Profiles
	
	DIR_CHECK_CREATE ${OUTPUTDIR}
	
	GENELIST_LIST=(up_up_genelist.txt down_down_genelist.txt)
	GENETYPE_LIST=(up_up down_down)

	
	for (( i = 0 ; i < ${#GENELIST_LIST[@]} ; i++ ))
	do
	GENELISTFILE=${GENELIST_LIST[$i]}
	GENETYPE=${GENETYPE_LIST[$i]}
	local GAP=400
	
	for CELLTYPE in WT_Treg_D20 Hdac12_KO_Treg_D20
	do
		for HISMODS in H3K9Ac H3K27Ac
		do

		if [ $CELLTYPE = WT_Treg ]
		then
		NORMALIZATION=1.0
		else
		NORMALIZATION=1.0
		fi

		READDIR=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Hdac/Treg_HistoneMarks_ChIPSeq_Sep2016/SICER/${CELLTYPE}_${HISMODS}
		READFILE=${CELLTYPE}_${HISMODS}-W200-G${GAP}-FDR0.0001-islandfiltered.bed
		OUTFILE=${OUTPUTDIR}/${CELLTYPE}_${HISMODS}_on_${GENETYPE}_${REGIONTYPE}.txt
		echo "python ${EXE_PATH} -b ${READDIR}/${READFILE} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} \
		-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTFILE}"
		python ${EXE_PATH} -b ${READDIR}/${READFILE} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} \
		-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTFILE}
		
		echo ""
		echo ""
		done
	done
done

	echo "done"


	}

#RUN_Two_Conditions_diagonal_plot(){}

##END OF FUNDEMENTAL FUNCTIONS

######################################

########################################################################

######################################

##	MODOLES
########################################################################
########################################################################

RUN_BOWTIE2(){
	### RUN_BOWTIE2 $1 $2
	local INPUT_NAME=${1}
	local SPECIES=${2}

	CHECK_arguments $# 2
	echo "RUN_BOWTIE2"
	local THREADS=4
	#### OUTPUT FORMAT
	local OUTPUT_BOWTIE2_FOLDER="${__EXE_PATH}/${INPUT_NAME}/bowtie2_results"
	DIR_CHECK_CREATE ${OUTPUT_BOWTIE2_FOLDER}
	local OUTPUT_BOWTIE2="${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.sam"
	
case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
	;;
	"mm10") 
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/HG19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
	;;
	"hg38")
	echo "Reference SPECIES is ${SPECIES}"
	local BOWTIEINDEXS="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/HG38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac

	#BOWTIEINDEXS_mm9_pMXs_combo="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/mm9_pMXs_combo/mm9_pMXs_combo"
	#BOWTIEINDEXS_mm9_combo_MMLV_pMXs="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/mm9_comdo_MMLV_pMXs/mm9_comdo_MMLV_pMXs"
	#BOWTIEINDEXS_MMLV_pMXs_combo="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/MMLV_pMXs_combo/MMLV_pMXs_combo"
	#BOWTIEINDEXS_MMLV="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/MMLV_index/MMLV"
	#BOWTIEINDEXS_mm10_pMXs_combo="/home/lxiang/Raw_Data/Paul/34bc/Bowtie2_Indexes/MM10_pMXs_combo/mm10_pMXs_combo"
	
########################################################################
	
	echo ""

	if [ -n "${__FASTQ_DIR_R1[0]}" -a -n "${__FASTQ_DIR_R2[0]}" ]
	then
	echo "Pair End Mode"
	echo "bowtie2 -p $THREADS -t --no-unal --non-deterministic -x $BOWTIEINDEXS -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S $OUTPUT_BOWTIE2 --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
	bowtie2 -p ${THREADS} -t --no-unal --non-deterministic -x ${BOWTIEINDEXS} -1 $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -2 $(echo ${__FASTQ_DIR_R2[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
	
	#### concordantly pair output
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x $BOWTIEINDEXS -1 ${__FASTQ_DIR_R1[0]} -2 ${__FASTQ_DIR_R2[0]} -S ${OUTPUT_BOWTIE2} --un-conc-gz ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R%.fastq.gz
	#### using un concordantly pair-ends do bowtie2 again.
	#echo "bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R1.fastq.gz -2 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R2.fastq.gz -S ${OUTPUT_BOWTIE2_unaligned}"
	#bowtie2 -p $THREADS --end-to-end --very-sensitive -k 1 --no-mixed --no-discordant --no-unal -x ${BOWTIEINDEXS} -1 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R1.fastq.gz -2 ${OUTPUT_BOWTIE2_FOLDER}/un_conc_aligned_R2.fastq.gz -S ${OUTPUT_BOWTIE2_unaligned}
	
	else
	echo "Single End Mode."
	echo "bowtie2 -p $THREADS --no-unal --non-deterministic -x $BOWTIEINDEXS -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S $OUTPUT_BOWTIE2"
	bowtie2 -p ${THREADS} --no-unal --non-deterministic -x ${BOWTIEINDEXS} -U $(echo ${__FASTQ_DIR_R1[*]} | tr " " ",") -S ${OUTPUT_BOWTIE2}
	fi

	echo ""
###sam2bam   ##FROM MACS2
	if [ ! -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
	echo "samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam"
	samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${OUTPUT_BOWTIE2} > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam 
	fi
### Then clear sam file.
	if [ -f ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ];then
	echo "rm ${OUTPUT_BOWTIE2}"
	rm ${OUTPUT_BOWTIE2}
	fi
###bam2bed

	echo "bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed"
	bamToBed -i ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam > ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bed

###bamToBedGraph

	#samtools sort ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}.bam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam


###sam2bed  Normal
#echo "samtools view $OUTPUT_BOWTIE2 -Sb | bamToBed -i stdin > $OUTPUT_BOWTIE2_FOLDER/$INPUT_NAME.bed"
#samtools view $OUTPUT_BOWTIE2 -Sb | bamToBed -i stdin > $OUTPUT_BOWTIE2_FOLDER/$INPUT_NAME.bed 


###bed2bedgraph
	#echo "bedtools genomecov -ibam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam -bga > ${INPUT_NAME}.bdg"
	#bedtools genomecov -ibam ${OUTPUT_BOWTIE2_FOLDER}/${INPUT_NAME}_sorted.bam -bga > ${INPUT_NAME}.bdg

	#RUN_BigGraph2BigWig ${OUTPUT_BOWTIE2_FOLDER} ${INPUT_NAME} "34bc"

	}

RUN_SICER(){
#### Usage: RUN_SICER $1 $2 ($1 is input for SICER. $2 is the CONTRO Library)
echo "RUN_SICER"
CHECK_arguments $# 2
local EXEDIR="/home/lxiang/Software/SICER1.1/SICER"

local GAP_SET=(400) # 600 400 200)
local REDUNDANCY=1
local INPUT_NAME=$1
local INPUI_CON=$2
local FRAGMENT_SIZE=1
local WINDOW_SIZE=200
local EFFECTIVEGENOME=0.85 
local FDR=0.05





local IN_SICER_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
local IN_SICER_FILES=${INPUT_NAME}.bed

local CONTRO_SICER_DIR=${__EXE_PATH}/${INPUI_CON}/bowtie2_results
local CONTRO_SICER_FILE=${INPUI_CON}.bed

local OUT_SICER_FOLDER=${__EXE_PATH}/${INPUT_NAME}/SICER_results
DIR_CHECK_CREATE ${OUT_SICER_FOLDER}

cd ${OUT_SICER_FOLDER}

echo "sh ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} mm9 ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log"
sh ${EXEDIR}/SICER.sh ${IN_SICER_FOLDER} ${IN_SICER_FILES} ${CONTRO_SICER_DIR} ${OUT_SICER_FOLDER} mm9 ${REDUNDANCY} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVEGENOME} ${GAP_SET} ${FDR} ${CONTRO_SICER_FILE} > ${OUT_SICER_FOLDER}/${INPUT_NAME}_SICER.log



RUN_Wig2BigWig ${OUT_SICER_FOLDER} ${INPUT_NAME}-W200-normalized
RUN_Wig2BigWig ${OUT_SICER_FOLDER} ${INPUI_CON}-W200-normalized


	}

RUN_MACS(){
#### kind outdated module
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

local INPUT_NAME=${1}
local INPUI_CON=${2}
local INPUT_LABEL="TCF1_032018"
#local INPUT_LABEL=${INPUT_NAME: 7:3}
##Skip out of Sample_ (7), and forward with 4 more digits.

#local FDR=0.05
local p_value=0.00001

local IN_FOLDER=${__EXE_PATH}/${INPUT_NAME}/bowtie2_results
local IN_FILES=${INPUT_NAME}.bam

local CONTRO_DIR=${__EXE_PATH}/${INPUI_CON}/bowtie2_results
local CONTRO_FILE=${INPUI_CON}.bam

local OUT_FOLDER=${__EXE_PATH}/${INPUT_NAME}/MACS2_results_BAMPE
DIR_CHECK_CREATE ${OUT_FOLDER}
cd ${OUT_FOLDER}

### --BAMPE
echo "python ${EXEDIR}/macs2 callpeak --format BAMPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value}" # -m 4  -q ${FDR}"
python ${EXEDIR}/macs2 callpeak --format BAMPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value} # -m 4 -q ${FDR}

### --BEDPE  !!!!! BEDPE is a new formate of pair end reads, do not treat bed format as bedpe!!!!!
#echo "python ${EXEDIR}/macs2 callpeak --format BEDPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value}" # -m 4  -q ${FDR}"
#python ${EXEDIR}/macs2 callpeak --format BEDPE -t ${IN_FOLDER}/${IN_FILES} -c ${CONTRO_DIR}/${CONTRO_FILE} --outdir ${OUT_FOLDER} -g 'mm' -n ${INPUT_NAME} -B -p ${p_value} # -m 4 -q ${FDR}
####      ${INPUT_NAME}_treat_pileup  for input

RUN_BigGraph2BigWig ${OUT_FOLDER} ${INPUT_NAME}_treat_pileup ${INPUT_LABEL}
RUN_BigGraph2BigWig ${OUT_FOLDER} ${INPUT_NAME}_control_lambda ${INPUT_LABEL}
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
local CON1_FILE=${CON1_TREAT}.bam

local CON1_CONTRO_FOLDER=${__EXE_PATH}/bowtie2_results
local CON1_CONTRO_FILE=${CON1_CONTRO}.bam

local CON2_FOLDER=${__EXE_PATH}/bowtie2_results
local CON2_FILE=${CON2_TREAT}.bam

local CON2_CONTRO_FOLDER=${__EXE_PATH}/bowtie2_results
local CON2_CONTRO_FILE=${CON2_CONTRO}.bam
	fi

########################################################################
########################################################################
local CON1_FOLDER=${__EXE_PATH}/${CON1_TREAT}/bowtie2_results
local CON1_FILE=${CON1_TREAT}.bam

local CON1_CONTRO_FOLDER=${__EXE_PATH}/${CON1_CONTRO}/bowtie2_results
local CON1_CONTRO_FILE=${CON1_CONTRO}.bam

local CON2_FOLDER=${__EXE_PATH}/${CON2_TREAT}/bowtie2_results
local CON2_FILE=${CON2_TREAT}.bam

local CON2_CONTRO_FOLDER=${__EXE_PATH}/${CON2_CONTRO}/bowtie2_results
local CON2_CONTRO_FILE=${CON2_CONTRO}.bam

########################################################################

########################################################################
local OUT_FOLDER=${__EXE_PATH}/MACS2_Diff_results/${CON1_NAME}_vs_${CON2_NAME}
DIR_CHECK_CREATE ${OUT_FOLDER}

echo "python ${EXEDIR}/macs2 callpeak --format BAMPE -B -t ${CON1_FOLDER}/${CON1_FILE} -c ${CON1_CONTRO_FOLDER}/${CON1_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON1_NAME} -g 'mm' --nomodel --extsize 200"
python ${EXEDIR}/macs2 callpeak --format BAMPE -B -t ${CON1_FOLDER}/${CON1_FILE} -c ${CON1_CONTRO_FOLDER}/${CON1_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON1_NAME} -g 'mm' -p ${p_value} # --nomodel --extsize 200 

echo "python ${EXEDIR}/macs2 callpeak --format BAMPE -B -t ${CON2_FOLDER}/${CON2_FILE} -c ${CON2_CONTRO_FOLDER}/${CON2_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON2_NAME} -g 'mm' --nomodel --extsize 200"
python ${EXEDIR}/macs2 callpeak --format BAMPE -B -t ${CON2_FOLDER}/${CON2_FILE} -c ${CON2_CONTRO_FOLDER}/${CON2_CONTRO_FILE} --outdir ${OUT_FOLDER} -n ${CON2_NAME} -g 'mm' -p ${p_value} # --nomodel --extsize 200 

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

RUN_CUFFDIFF(){
	### RUN_CUFFDIFF $1 $2
	####
CHECK_arguments $# 12
local INPUT_Args=("$@")
local SAMPLE_NUM=${#INPUT_Args[*]}
local GTFFILE=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2014-05-23-16-05-24/Genes/genes.gtf
local OUTPUT_Cuffdiff=${__EXE_PATH}/Cuffdiff_Results
local THREADS=4
DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}

#### Preparing the .bam files for cuffdiff
for (( i = 0; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do	
	local DATA_BAM[$i]=${__EXE_PATH}/${INPUT_Args[i]}/Tophat_results/${INPUT_Args[i]}.bam
	echo ${DATA_BAM[$i]}
done

### Four GROUP, each group has 3 DATA FILES, 
########################################################################
if [ ${SAMPLE_NUM} -eq "12" ]
then

	DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}/A_vs_B
	echo "${INPUT_Args} Data files are loading..."
	echo "cuffdiff -o ${OUTPUT_Cuffdiff}/A_vs_B -p ${THREADS} -L "Group_A","Group_B" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]},${DATA_BAM[2]} ${DATA_BAM[3]},${DATA_BAM[4]},${DATA_BAM[5]}"
	cuffdiff -o ${OUTPUT_Cuffdiff}/A_vs_B -p ${THREADS} -L "Group_A","Group_B" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]},${DATA_BAM[2]} ${DATA_BAM[3]},${DATA_BAM[4]},${DATA_BAM[5]}


	DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}/C_vs_D
	echo "cuffdiff -o${OUTPUT_Cuffdiff}/C_vs_D -p ${THREADS} -L "Group_C","Group_D" ${GTFFILE} ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]} ${DATA_BAM[9]},${DATA_BAM[10]},${DATA_BAM[11]}"
	cuffdiff -o${OUTPUT_Cuffdiff}/C_vs_D -p ${THREADS} -L "Group_C","Group_D" ${GTFFILE} ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]} ${DATA_BAM[9]},${DATA_BAM[10]},${DATA_BAM[11]}

	DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}/A_vs_C
	echo "cuffdiff -o ${OUTPUT_Cuffdiff}/A_vs_C -p ${THREADS} -L "Group_A","Group_C" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]},${DATA_BAM[2]} ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]}"
	cuffdiff -o ${OUTPUT_Cuffdiff}/A_vs_C -p ${THREADS} -L "Group_A","Group_C" ${GTFFILE} ${DATA_BAM[0]},${DATA_BAM[1]},${DATA_BAM[2]} ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]}


	DIR_CHECK_CREATE ${OUTPUT_Cuffdiff}/B_vs_D
	echo "cuffdiff -o ${OUTPUT_Cuffdiff}/B_vs_D -p ${THREADS} -L "Group_B","Group_D" ${GTFFILE}  ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]} ${DATA_BAM[9]},${DATA_BAM[10]},${DATA_BAM[11]}"
	cuffdiff -o ${OUTPUT_Cuffdiff}/B_vs_D -p ${THREADS} -L "Group_B","Group_D" ${GTFFILE}  ${DATA_BAM[6]},${DATA_BAM[7]},${DATA_BAM[8]} ${DATA_BAM[9]},${DATA_BAM[10]},${DATA_BAM[11]}

fi
########################################################################

echo "CuffDiff Completed!"
	}

RUN_TOPHAT (){
#### Usage: RUN_TOPHAT $1 $2 $3 $4

echo "RUN_TOPHAT_ANALYSIS"
CHECK_arguments $# 4

### Operation PARAMETERS Setting
local INPUT_NAME=${1}
local PROJECT_NAME=${2}
local SPECIES=${3}
local Data_Provider=${4}

#### OUTPUT FORMAT
local OUTPUT_TOPHAT_FOLDER="${__EXE_PATH}/${INPUT_NAME}/Tophat_results"
DIR_CHECK_CREATE ${OUTPUT_TOPHAT_FOLDER}

local WINDOW_SIZE=200
local FRAGMENT_SIZE=50
local THREADS=4
########################################################################
case ${SPECIES} in
	"mm9")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2014-05-23-16-05-24/Genes/genes.gtf
	local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
	;;
	"hg19")
	echo "Reference SPECIES is ${SPECIES}"
	local GTFFILE=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf
	local BOWTIEINDEXS=/home/data/Annotation/iGenomes/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	;;
	*)
	echo "ERR: Did Not Find Any Matched Reference...... Exit"
	exit
	;;
esac
########################################################################
local SICER_DIR=/home/lxiang/Software/SICER1.1/SICER
local EXEDIR=$SICER_DIR/extra/tools/wiggle
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
	
	echo "bamToBed -i ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam > ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed"
	bamToBed -i ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bam > ${OUTPUT_TOPHAT_FOLDER}/${INPUT_NAME}.bed
	
	echo "sh $EXEDIR/bed2wig.sh ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}"
	sh $EXEDIR/bed2wig.sh ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME} ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${SPECIES}
	
	echo "RUN_Wig2BigWig ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME}-W${WINDOW_SIZE}-RPKM ${PROJECT_NAME} ${SPECIES} ${Data_Provider}"
	RUN_Wig2BigWig ${OUTPUT_TOPHAT_FOLDER} ${INPUT_NAME}-W${WINDOW_SIZE}-RPKM ${PROJECT_NAME} ${SPECIES} ${Data_Provider}
	
	echo "One TopHat is completed."
	echo ""
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
	CHECK_arguments $# 1
	
	EXEDIR=/home/lxiang/Software/python_tools

	GTFDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gtf_files
	GTFFILE=mm9_genes.gtf
### This parameter means that Promoter [-1k TSS, +1k TSS] Gene_body[+1k TSS. TES]
	PROMOTER_UPSTREAM_EXTENSION=1001   # Upstream extension, the distance from TSS.
	TSS_REGION_LENGTH=2000


	INPUTFILE=${1}.bed
	OUTPUTFILE=${1}_distribution.txt

	echo "python $EXEDIR/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE"
	python $EXEDIR/peaks_count_on_genic_region.py -i ${__RAW_DATA_PATH_DIR}/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o ${__EXE_PATH}/$OUTPUTFILE
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

EMAIL_ME(){
	CHECK_arguments $# 1
	echo "Start at ${1} " | mailx -v -s "Project: + ${2} Finished" lux@gwu.edu
	}

FUNC_Download (){
	CHECK_arguments $# 1
### Download Login information and the download directory.
#### -nH --cut-dirs=3   Skip 3 directory components.
	local Web_Address="https://jumpgate.caltech.edu/runfolders/volvox02/180426_SN787_0801_BHGW2YBCX2/Unaligned"
	local Directory_Skip_Num=4
	local USER="gec"
	local PASSWORD="aardvark dryer rummage"
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
	echo "Input Dir: $@"
	local Dirs=$@
	for solo_dir in $Dirs
	do
		if [ ! -d $solo_dir ];then
			echo "Dir check and create is $solo_dir"
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
