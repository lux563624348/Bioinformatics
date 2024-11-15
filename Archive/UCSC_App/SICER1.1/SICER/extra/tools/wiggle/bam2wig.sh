#!/bin/bash

SICER=~/Software/SICER1.1/SICER

SAMPLEDIR=$1
SAMPLE=$2
WINDOW_SIZE=$3
FRAGMENT_SIZE=$4
SPECIES=$5

EXEDIR=$SICER/extra/tools/wiggle

cd $SAMPLEDIR

#bamToBed is from BEDTools
echo "bamToBed -i $SAMPLE.bam > $SAMPLE.bed"
echo ""
bamToBed -i $SAMPLE.bam > $SAMPLE.bed
echo "sh $EXEDIR/bed2wig.sh $SAMPLEDIR $SAMPLE $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES"
echo ""
sh $EXEDIR/bed2wig.sh $SAMPLEDIR $SAMPLE $WINDOW_SIZE $FRAGMENT_SIZE $SPECIES

