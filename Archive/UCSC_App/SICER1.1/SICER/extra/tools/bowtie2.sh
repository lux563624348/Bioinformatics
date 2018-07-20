#\!/bin/bash

#Run Bowtie2 on paired-end RNA Seq data
#Generate BED format for downstream analysis, need samtools and BedTools
#Generate wig files for browser
#zip the huge files

SEQ1=${1}_1.fq
SEQ2=${1}_2.fq
OUTPUT=$1.sam

## indexes can be downloaded from bowtie website and I put them in the following folder
#BOWTIEINDEXES="/home/wpeng/Downloads/bowtie_annotation"
BOWTIEINDEXES=$2
#SPECIES=mm9
SPECIES=$3

BOWTIE2="/home/wpeng/Downloads/bowtie2-2.0.0-beta5/bowtie2"
SEQDIR=.

# trim 1 bp from 5 prime because of the N there in lane 2 
echo "$BOWTIE2 -5 1 -p 4 -t -x $BOWTIEINDEXES/$SPECIES -1 $SEQDIR/$SEQ1 -2 $SEQDIR/$SEQ2 -S $OUTPUT"
$BOWTIE2 -5 1 -p 4 -t -x $BOWTIEINDEXES/$SPECIES  -1 $SEQDIR/$SEQ1 -2 $SEQDIR/$SEQ2 -S $OUTPUT

#sam2bed
samtools view $OUTPUT -Sb | bamToBed -i stdin > $1.bed 

#bed2wigs
EXEDIR=/home/data/SICER1.1/SICER/extra/tools/wiggle
SAMPLEDIR=.
SAMPLE=$1
WINDOW_SIZE=10
FRAGMENT_SIZE=0
sh $EXEDIR/bed2wig $SAMPLEDIR $SAMPLE $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES

#zip original fq files
gzip $SEQ1
gzip $SEQ2
gzip $OUTPUT
