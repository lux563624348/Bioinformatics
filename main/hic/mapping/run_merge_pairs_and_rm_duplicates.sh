#!/bin/bash

GENOME_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

INPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905/iterative_mapping
OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905/iterative_mapping

for SAMPLE in WT_na_CD8 Tcf1_KO_na_CD8
do

	echo "python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE"
	python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE

done

echo "done"
