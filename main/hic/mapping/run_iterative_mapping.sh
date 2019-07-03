#!/bin/bash

EXEDIR=/home/lxiang/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/HiC/1905

BOWTIE_PATH=/usr/local/bowtie2-2.2.6/bin/bowtie2
BOWTIE_INDEX_PATH=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/UCSC/Mouse_Genome/MM9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

OUTPUT_BAM_DIR=/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/CD8-HP/HiC/1905

MIN_SEQ_LEN=25
LEN_STEP=5

SEQ_START=0
SEQ_END=50

FASTQ_DIR=/home/lxiang/cloud_research/PengGroup/XLi/Raw_Data/Haihui/CD8-HP/HiC/1905

for PAIR in {1..2}
do
	for CELL_TYPE in WT_na_CD8 Tcf1_KO_CD8
	do

		if [ $CELL_TYPE = WT_na_CD8 ]
		then
			FILENAME_SET=$(cd $FASTQ_DIR; find . -name "WT_na_CD8*R${PAIR}*" | xargs)	
		else
			FILENAME_SET=$(cd $FASTQ_DIR; find . -name "Tcf1_KO_CD8*R${PAIR}*" | xargs)
		fi

		IFS=' ' read -a FILENAME_LIST <<< "$FILENAME_SET"
		
		SAMPLE=${CELL_TYPE}_R${PAIR}

		FASTQ_PATH=$FASTQ_DIR/${FILENAME_LIST[$i]}
		OUTPUT_BAM_PATH=$OUTPUT_BAM_DIR/$CELL_TYPE/${SAMPLE}.bam

		echo "python $EXEDIR/run_iterative_mapping.py -b $BOWTIE_PATH -i $BOWTIE_INDEX_PATH -q $FASTQ_PATH -o $OUTPUT_BAM_PATH -m $MIN_SEQ_LEN -t $LEN_STEP -s $SEQ_START -e $SEQ_END"
		python $EXEDIR/run_iterative_mapping.py -b $BOWTIE_PATH -i $BOWTIE_INDEX_PATH -q $FASTQ_PATH -o $OUTPUT_BAM_PATH -m $MIN_SEQ_LEN -t $LEN_STEP -s $SEQ_START -e $SEQ_END
		echo ""
		echo ""

	done
done

echo "done"
