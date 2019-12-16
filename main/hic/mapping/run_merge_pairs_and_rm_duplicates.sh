#!/bin/bash

GENOME_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

INPUTDIR=/home/lxiang/Data_Processing/Haihui/CD8-HP/HiC_201910/Naive
OUTPUTDIR=/home/lxiang/Data_Processing/Haihui/CD8-HP/HiC_201910/Naive

__INPUT_SAMPLE_SET=(
nWT3_Rep1
nWT3_Rep2
nDKO2_Rep1
nDKO2_Rep2
)

for SAMPLE in ${__INPUT_SAMPLE_SET[*]}
do

	echo "python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE"
	python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE

done

echo "done"
