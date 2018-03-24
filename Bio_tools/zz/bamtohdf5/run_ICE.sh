#!/bin/bash

GENOME_PATH=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping
RESOLUTION=10000

for SAMPLE in stim_WT_CD8
do

	echo "python run_ICE.py -g $GENOME_PATH -n $SAMPLE -o $OUTPUTDIR -r $RESOLUTION"
	python run_ICE.py -g $GENOME_PATH -n $SAMPLE -o $OUTPUTDIR -r $RESOLUTION
	echo ""
	echo ""

done

echo "done"
