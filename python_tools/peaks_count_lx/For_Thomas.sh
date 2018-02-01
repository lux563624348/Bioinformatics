#!/bin/bash

#	nohup bash run_main.sh 1>default.log 2>default_err.log &
########################################################################
## 9/29/2017
## By Xiang Li, basing on Zhouhao's scripts.
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


#			NEED MODIFICATION FOR DIFFERENT PROJECT

########################################################################
## GLOBAL VARIABLES
########################################################################
__RAW_DATA_PATH_DIR=$(pwd)
#### Execution or Output directory
__EXE_PATH=$(pwd)/OUTPUT
########################################################################

__INPUT_SAMPLE_DIR_List=(
gene_promoterGenebody_U_4k_iv_unique
INPUT_A
)



echo "INPUT_SAMPLE_DIR_List= (${__INPUT_SAMPLE_DIR_List[*]})"
SAMPLE_NUM=${#__INPUT_SAMPLE_DIR_List[*]}
########################################################################

########################################################################
##	MAIN BODY
########################################################################

main() {
#### Saving DIR Check and Create
DIR_CHECK_CREATE $__RAW_DATA_PATH_DIR
DIR_CHECK_CREATE $__EXE_PATH



echo ""

for (( i = 1; i <= $(expr $SAMPLE_NUM - 1); i++ ))
do
	
	echo "RUN_bed_intersect ${__INPUT_SAMPLE_DIR_List[1]} ${__INPUT_SAMPLE_DIR_List[0]}"

	RUN_bed_intersect ${__INPUT_SAMPLE_DIR_List[0]} ${__INPUT_SAMPLE_DIR_List[1]}
	
	cd $__EXE_PATH
	rm genes.bed peaks.bed peaks_counting.bed 
	break
done



	echo "End Date: `date`"

}

########################################################################
##	MODULES
########################################################################

RUN_bed_intersect(){
	#### usage RUN_bed_intersect $1 $2

	local Input_A1=${__RAW_DATA_PATH_DIR}/${1}.bed
	local Input_B1=${__RAW_DATA_PATH_DIR}/${2}.bed
	cd ${__EXE_PATH}
	local Title=""
	
	####################################################################
	echo "$Title"
	echo ""
	
	
	echo "Input A1"
	wc -l ${Input_A1}
	echo "Input B1"
	wc -l ${Input_B1}
	


	Output_Path="${__EXE_PATH}"
	
	bedtools sort -i ${__RAW_DATA_PATH_DIR}/${2}.bed > ${__RAW_DATA_PATH_DIR}/${2}_sorted.bed
	
	Input_B1=${__RAW_DATA_PATH_DIR}/${2}_sorted.bed
	
	echo "bedtools intersect -wa -a ${Input_A1} -b ${Input_B1} -c > ${Output_Path}/genes.bed"
	bedtools intersect -wa -a ${Input_A1} -b ${Input_B1} -c > ${Output_Path}/genes.bed
	#wc -l ${Output_Path}/genes.bed
######This line is specific for NANPING 10.22.17

	echo "bedtools intersect -u -a ${Input_B1} -b ${Input_A1} > ${Output_Path}/peaks.bed"
	bedtools intersect -u -a ${Input_B1} -b ${Input_A1} > ${Output_Path}/peaks.bed
	
	echo "bedtools intersect -a ${Input_B1} -b ${Input_A1} -c > ${Output_Path}/peaks_counting.bed"
	bedtools intersect -a ${Output_Path}/peaks.bed -b ${Input_A1} -c > ${Output_Path}/peaks_counting.bed
	#wc -l ${Output_Path}/peaks.bed
#####This line is specific for NANPING 10.22.17

	python ${__RAW_DATA_PATH_DIR}/read_count.py ${__RAW_DATA_PATH_DIR} ${__EXE_PATH} ${2}
	
	
	
	}

DIR_CHECK_CREATE (){
### Saving DIR Check and Create
	echo "Dir check and create are: $@"
	if [ ! -d $@ ];then
	mkdir -p $@
	fi
	}


##### Following Line is very IMPORTANT 
main "$@"

