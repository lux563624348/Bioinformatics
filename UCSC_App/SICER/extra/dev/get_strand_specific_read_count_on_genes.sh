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


#UCSCDIR=/home/data/hg18/CD16/
#UCSCFILE=hg18-Agilent-19554-genes.ucsc
UCSCDIR=/home/wyang/data/hg18/annotation/
UCSCFILE=knownGene.txt

# LIBA=72h_EM_K27
# READDIR=/home/data/hg18/CD8/2010data/processed
# READFILE=$LIBA-W200-G600-E500-islandfiltered.bed 

LIBA=GA2068-hg18-CD16-A1-RNAseq
READDIR=/home/data/hg18/CD16/raw/set1
READFILE=$LIBA.bed

FRAGMENTSIZE=0
#TYPE=GeneBody
TYPE=ExonicRegion
UPSTREAM=0
DOWNSTREAM=0
OUTFILE=$LIBA-Exon-ReadCount.dat

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=get_strand_specific_read_count_on_genes.py
OUTFILE=${LIBA}-${TYPE}-ReadCount.dat
echo "python $EXEDIR/$EXE   -b $READDIR/$READFILE -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE"
python $EXEDIR/$EXE   -b $READDIR/$READFILE -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE



#echo "python $EXEDIR/$EXE   -b $READDIR/$READFILE -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE"
#python $EXEDIR/$EXE   -b $READDIR/$READFILE -f $FRAGMENTSIZE -g $UCSCDIR/$UCSCFILE -r $TYPE -u $UPSTREAM -d $DOWNSTREAM -o $OUTFILE