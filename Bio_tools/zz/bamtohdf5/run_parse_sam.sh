#!/bin/bash

GENOME_PATH=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping

for SAMPLE in stim_WT_CD8 stim_DKO_CD8
do

	INPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/${SAMPLE}

	for PART in p1 p2 p3 p4
	do

		INPUTFILE_R1=${SAMPLE}_R1_${PART}.bam
		INPUTFILE_R2=${SAMPLE}_R2_${PART}.bam
 
		OUTPUTFILE=${SAMPLE}_${PART}_mapped_reads.hdf5

		echo "python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE"
		python run_parse_sam.py -g $GENOME_PATH -1 $INPUTDIR/$INPUTFILE_R1 -2 $INPUTDIR/$INPUTFILE_R2 -o $OUTPUTDIR/$OUTPUTFILE 
		echo ""
		echo ""
		
	done
done

echo "done"
