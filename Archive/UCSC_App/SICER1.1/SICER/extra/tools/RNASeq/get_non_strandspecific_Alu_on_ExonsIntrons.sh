#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq
ALUDIR=/home/data/hg19/Annotation/RepElements
ALUS=AluElements_hg19.txt
#ALUS=AluElements_hg19_sample10000.txt
GENESDIR=$EXEDIR
GENES=hg19_EntrezID_filtered_samestrand_collisonremoved_simplegene.pkl
SPECIES=hg19

echo "python $EXEDIR/get_non_strandspecfic_Alu_on_ExonsIntrons.py -a $ALUDIR/$ALUS -u $GENESDIR/$GENES  -s $SPECIES -o simplegenes_$SPECIES"
python $EXEDIR/get_non_strandspecfic_Alu_on_ExonsIntrons.py -a $ALUDIR/$ALUS -u $GENESDIR/$GENES  -s $SPECIES -o simplegenes_$SPECIES