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


#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq
EXE=get_strandspecific_read_count_on_genes.py

READSDIR=/home/data/hg19/JunZhu/raw/CD4/Round2

ANNOTATIONDIR=/home/data/hg19/Annotation
ANNOTATION=refFlat_hg19_EntrezID_filtered.ucsc

TYPE=GeneBody
UPSTREAM=0
DOWNSTREAM=1000
FRAGMENTSIZE=0

for READS in Rest_R1_unique Active_R1_unique
do
READFILEF=${READS}_P.bed
READFILER=${READS}_N.bed

OUTFILE=$READS-on-$TYPE.dat

echo "python $EXEDIR/get_strandspecific_read_count_on_genes.py -p $READSDIR/$READFILEF -n $READSDIR/$READFILER -f $FRAGMENTSIZE -g $ANNOTATIONDIR/$ANNOTATION  -r $TYPE -u $UPSTREAM -d $DOWNSTREAM  -o $OUTFILE -s hg19 "

python $EXEDIR/get_strandspecific_read_count_on_genes.py -p $READSDIR/$READFILEF -n $READSDIR/$READFILER -f $FRAGMENTSIZE -g $ANNOTATIONDIR/$ANNOTATION  -r $TYPE -u $UPSTREAM -d $DOWNSTREAM  -o $OUTFILE -s hg19 

done

