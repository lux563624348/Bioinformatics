#!/bin/bash
# Authors: Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).


##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

#LIBA=72h_EM_K27
#READDIR=/home/data/hg18/CD8/2010data/processed
#READFILE=$LIBA-W200-G600-E500-islandfiltered.bed 


ISLANDDIR=/home/data/SICER1.1/SICER/extra/test
#ISLANDFILE=$LIBA-W200-G600-E500.scoreisland
ISLANDFILE=testisland

READDIR=$1
LIBA=$2


#LIBA=H3K27me3-1-removed
#READDIR=/home/data/hg18/CD4/processed 

#LIBA=H3K27me3
#READDIR=/home/data/hg18/CD4 
READFILE=$LIBA.bed


EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=calculate_correlation_dev.py


TYPE=cross
SPECIES=hg18
BINSIZE=10
MINPOINTS=5
MAXDISTANCE=100
RESOLUTION=1
SHIFT=0
OUTFILE=${LIBA}-${TYPE}-B$BINSIZE-R$RESOLUTION-S$SHIFT-correlation

echo "python $EXEDIR/$EXE  -b $READDIR/$READFILE  -s $SPECIES -i $ISLANDDIR/$ISLANDFILE -w $BINSIZE -m $MINPOINTS -n $MAXDISTANCE -r $RESOLUTION -t $TYPE -f $SHIFT -o $OUTFILE"


python $EXEDIR/$EXE  -b $READDIR/$READFILE  -s $SPECIES -i $ISLANDDIR/$ISLANDFILE -w $BINSIZE -m $MINPOINTS -n $MAXDISTANCE -r $RESOLUTION -t $TYPE -f $SHIFT -o $OUTFILE

