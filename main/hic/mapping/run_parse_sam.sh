#!/bin/bash

GENOME_PATH=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

OUTPUTDIR=/home/lxiang/Data_Processing/Haihui/CD8-HP/HiC_201910/Naive

__INPUT_SAMPLE_SET=(
19092FL-06-01_S1_L001
19092FL-06-01_S1_L002
19092FL-06-02_S2_L003
19092FL-06-02_S2_L004
)

for SAMPLE in ${__INPUT_SAMPLE_SET[*]}
do

	INPUTDIR=/home/lxiang/Data_Processing/Haihui/CD8-HP/HiC_201910/Naive/${SAMPLE}

	INPUTFILE_R1=${SAMPLE}_R1.bam
	INPUTFILE_R2=${SAMPLE}_R2.bam

	OUTPUTFILE=${SAMPLE}_mapped_reads.hdf5

	echo "python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE"
	python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE 
	echo ""
	echo ""
		
done

echo "done"
