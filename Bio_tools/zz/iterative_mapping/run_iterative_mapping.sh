#!/bin/bash

EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017

BOWTIE_PATH=/usr/local/bowtie2-2.2.6/bin/bowtie2
BOWTIE_INDEX_PATH=/home/data/Annotation/iGenomes/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

OUTPUT_BAM_DIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping

MIN_SEQ_LEN=25
LEN_STEP=5

SEQ_START=0
SEQ_END=50

FASTQ_DIR=/home/zzeng/cloud_research/PengGroup/ZZeng/raw_data/Haihui/Tcf1/Tcf1_HiCSeq_Jul2017

for PAIR in {1..2}
do
	for CELL_TYPE in stim_WT_CD8 stim_DKO_CD8
	do

		if [ $CELL_TYPE = stim_WT_CD8 ]
		then
			FILENAME_SET=$(cd $FASTQ_DIR; find . -name "WtCD8stim*R${PAIR}*" | xargs)	
		else
			FILENAME_SET=$(cd $FASTQ_DIR; find . -name "dKOCD8stm*R${PAIR}*" | xargs)
		fi

		IFS=' ' read -a FILENAME_LIST <<< "$FILENAME_SET"
		
		for i in {0..3}
		do
			SAMPLE=${CELL_TYPE}_R${PAIR}_p$(($i+1))
			
			FASTQ_PATH=$FASTQ_DIR/${FILENAME_LIST[$i]}
			OUTPUT_BAM_PATH=$OUTPUT_BAM_DIR/$CELL_TYPE/${SAMPLE}.bam

			echo "python $EXEDIR/run_iterative_mapping.py -b $BOWTIE_PATH -i $BOWTIE_INDEX_PATH -q $FASTQ_PATH -o $OUTPUT_BAM_PATH -m $MIN_SEQ_LEN -t $LEN_STEP -s $SEQ_START -e $SEQ_END"
			python $EXEDIR/run_iterative_mapping.py -b $BOWTIE_PATH -i $BOWTIE_INDEX_PATH -q $FASTQ_PATH -o $OUTPUT_BAM_PATH -m $MIN_SEQ_LEN -t $LEN_STEP -s $SEQ_START -e $SEQ_END
			echo ""
			echo ""

		done

	done
done

echo "done"
