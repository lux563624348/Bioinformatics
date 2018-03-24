#!/bin/bash

GENOME_PATH=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

INPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping
OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping

for SAMPLE in stim_WT_CD8
do

	echo "python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE"
	python run_merge_pairs_and_rm_duplicates.py -g $GENOME_PATH -i $INPUTDIR -o $OUTPUTDIR -n $SAMPLE

done

echo "done"
