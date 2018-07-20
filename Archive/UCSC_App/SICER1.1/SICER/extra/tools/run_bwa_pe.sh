#\!/bin/bash

#building index
#only need to run once
#bwa index -p mm9bwaidx -a bwtsw wg.fa

#alignment:
GENOMEDIR=/home/data/mm9
GENOME=mm9bwaidx

RAWDIR=/home/data/mm9/Lin/raw/Jan2013/Sample_LHPL001A
RAW=LHPL001A_ATCACG_L002_001
RAW1=LHPL001A_ATCACG_L002_R1_001
RAW2=LHPL001A_ATCACG_L002_R2_001

RAWFILE1=$RAW1.fastq.gz
RAWFILE2=$RAW2.fastq.gz

# -t using multiple threads
bwa aln -t 2  $GENOMEDIR/$GENOME $RAWDIR/$RAWFILE1 > $RAW1.sai
bwa aln -t 2  $GENOMEDIR/$GENOME $RAWDIR/$RAWFILE2 > $RAW2.sai 

#mapping (single-end use "samese"; pair-end use "sampe"
#--samse
#bwa samse /path/to/reference/ref.fa file1.sai /path/to/fastq/file1.txt  > file1.sam
#--sampe
bwa sampe $GENOMEDIR/$GENOME $RAW1.sai $RAW2.sai $RAWDIR/$RAWFILE1 $RAWDIR/$RAWFILE2 > $RAW.sam

#sam2bam (only keep uniquely mappable reads: -q 1)
UNIQUEREADS=${RAW}_unique_reads
MULTIPLEREADS=${RAW}_multiple_reads
samtools view -bS -q 1 $RAW.sam > $UNIQUEREADS.bam
samtools view -bS $RAW.sam > $MULTIPLEREADS.bam

#bam2bed
bamToBed -i $UNIQUEREADS.bam > $UNIQUEREADS.bed
bamToBed -i $MULTIPLEREADS.bam > $MULTIPLEREADS.bed