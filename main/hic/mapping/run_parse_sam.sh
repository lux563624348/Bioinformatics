#!/bin/bash

GENOME_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905/iterative_mapping

for SAMPLE in WT_na_CD8 Tcf1_KO_na_CD8
do

	INPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905/iterative_mapping/${SAMPLE}

	INPUTFILE_R1=${SAMPLE}_R1.bam
	INPUTFILE_R2=${SAMPLE}_R2.bam

	OUTPUTFILE=${SAMPLE}_mapped_reads.hdf5

	echo "python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE"
	python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE 
	echo ""
	echo ""
		
done

echo "done"
