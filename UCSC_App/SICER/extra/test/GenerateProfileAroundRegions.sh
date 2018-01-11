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

LIBA=H3K36ac_hg19
UCSCDIR=/home/data/hg19/Annotation
UCSCFILE=refFlat_hg19_sorted.ucsc

READDIR=.
READFILE=$LIBA-W200-G600-FDR.01-islandfiltered.bed

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=GenerateProfileAroundRegions.py
#EXE=plottingtest.py

TYPE=GENE
SPECIES=hg19
UPSTREAM=5000
DOWNSTREAM=5000
RESOLUTION=1000
WINDOWSIZE=1000
GENICPARTITION=40
PLUSSHIFT=72
MINUSSHIFT=72
OUTFILE=${LIBA}_on_${TYPE}_profile

echo "python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $READDIR/$READFILE -c $TYPE -o $OUTFILE -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -g $GENICPARTITION -p $PLUSSHIFT -m $MINUSSHIFT"

python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $READDIR/$READFILE -c $TYPE -o $OUTFILE -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -g $GENICPARTITION -p $PLUSSHIFT -m $MINUSSHIFT


