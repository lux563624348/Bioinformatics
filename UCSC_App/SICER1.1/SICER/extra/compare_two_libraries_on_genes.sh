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

LIBADIR=/home/data/hg18/CD8/2010data/raw
LIBA=0N_EZh2.bed

LIBBDIR=/home/data/hg18/CD8/2010data/raw
LIBB=0NEZh2.bed

UCSCDIR=/home/data/hg18/CD8
UCSCFILE=hg18-Agilent-19554-genes.ucsc

READDIR=/home/data/hg18/CD8/2010data/processed
READFILE=$LIBA-W200-G600-E500-islandfiltered.bed 

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=compare_two_libraries_on_genes.py

FRAGMENTSIZE=150
TYPE=Promoter
UPSTREAM=5000
DOWNSTREAM=5000
OUTFILE=${LIBA}-vs-$LIBB-on_${TYPE}.txt

echo "python $EXEDIR/$EXE  -a $LIBADIR/$LIBA -b $LIBBDIR/$LIBB -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE"

python $EXEDIR/$EXE  -a $LIBADIR/$LIBA -b $LIBBDIR/$LIBB -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE

