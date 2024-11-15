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

LIBA=72h_EM_K27
ISLANDDIR=/home/data/SICER1.1/SICER/extra/test
ISLANDFILE=testisland

READDIR=/home/data/hg18/CD8/2010data/processed
READFILE=$LIBA-W200-G600-E500-islandfiltered.bed 

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=filter_raw_tags_by_islands_dev.py

SHIFT=75
SPECIES=hg18
TYPE=Promoter
UPSTREAM=5000
DOWNSTREAM=5000
OUTFILE=${LIBA}-island-filtered-reads.bed

echo "python $EXEDIR/$EXE  -s $SPECIES -a $READDIR/$READFILE -i $SHIFT  -b $ISLANDDIR/$ISLANDFILE -o $OUTFILE"

python $EXEDIR/$EXE  -s $SPECIES -a $READDIR/$READFILE -i $SHIFT  -b $ISLANDDIR/$ISLANDFILE -o $OUTFILE
