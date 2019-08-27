#!/bin/bash

GENOME_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905/iterative_mapping/
RESOLUTION=10000

for SAMPLE in WT_na_CD8 Tcf1_KO_na_CD8
do

	echo "python run_ICE.py -g $GENOME_PATH -n $SAMPLE -o $OUTPUTDIR -r $RESOLUTION"
	python run_ICE.py -g $GENOME_PATH -n $SAMPLE -o $OUTPUTDIR -r $RESOLUTION
	echo ""
	echo ""

done

echo "done"
