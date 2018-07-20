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

N=100

LIBA=72h_EM_K27
#UCSCDIR=/home/data/hg18/CD8
#UCSCFILE=hg18-Agilent-19554-genes.ucsc
UCSCDIR=.
UCSCFILE=hg18-Agilent-19554-genes_head$N.ucsc

READDIR=/home/data/hg18/CD8/2010data/processed
READFILE=$LIBA-W200-G600-E500-islandfiltered.bed 

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=GenerateProfileMatrixAroundLocations.py

TYPE=TSS
SPECIES=hg18
UPSTREAM=5000
DOWNSTREAM=5000
RESOLUTION=100
WINDOWSIZE=500
NORMALIZATION=1
PLUSSHIFT=0
MINUSSHIFT=0
OUTFILE=${LIBA}_on_${TYPE}_R${RESOLUTION}_W${WINDOWSIZE}_head${N}_profileMatrix

echo "python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $READDIR/$READFILE -c $TYPE -o $OUTFILE  -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -p $PLUSSHIFT -m $MINUSSHIFT"

python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $READDIR/$READFILE -c $TYPE -o $OUTFILE  -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -p $PLUSSHIFT -m $MINUSSHIFT

