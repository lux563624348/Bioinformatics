#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq

#READSDIR=/home/data/hg19/JunZhu/raw/CD4/Round2
READSDIR = ./test	

for READS in Rest_R1_n100000_unique 
do
READFILEF=${READS}_P.bed
READFILER=${READS}_N.bed
FRAGMENTSIZE=0

ANNOTATIONDIR=/home/data/hg19/Annotation
ANNOTATION=hg19_EntrezID_filtered_samestrand_collisonremoved.pkl 

OUTFILE=$READS-on-EligibleEntrezGenes.dat

echo "python $EXEDIR/get_strandspecific_read_count_on_ExonsIntrons.py -f $READSDIR/$READFILEF -r $READSDIR/$READFILER -g $FRAGMENTSIZE -u $ANNOTATIONDIR/$ANNOTATION -o $OUTFILE -s hg19 "
python $EXEDIR/get_read_count_on_ExonsIntrons.py -f $READSDIR/$READFILEF -r $READSDIR/$READFILER -g $FRAGMENTSIZE -u $ANNOTATIONDIR/$ANNOTATION -o $OUTFILE -s hg19


# echo "sort -g -k1 $OUTFILE > $READS-on-EligibleEntrezGenes_sorted.dat"
# sort -g -k1 $OUTFILE > $READS-on-EligibleEntrezGenes_sorted.dat
# echo "mv $READS-on-EligibleEntrezGenes_sorted.dat > $OUTFILE"
# mv $READS-on-EligibleEntrezGenes_sorted.dat > $OUTFILE

echo ""
echo ""

done

