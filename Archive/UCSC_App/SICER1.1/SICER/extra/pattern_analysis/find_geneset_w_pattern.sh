#!/bin/bash
SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib:$SICER/extra
export PYTHONPATH


EXEDIR=/home/data/SICER1.1/SICER/extra/pattern_analysis
COMPLETEPATTERNFILE=PromoterDigitalPattern.txt
ANNOTATIONFILE=/home/data/mm9/Annotation/mm9_rmsk.ucsc
PATTERN=$1
OUTPUTFILE=${PATTERN}.ucsc

python $EXEDIR/find_gene_set_with_pattern.py -i $COMPLETEPATTERNFILE -k $ANNOTATIONFILE -p $PATTERN -o $OUTPUTFILE
echo "python $EXEDIR/find_gene_set_with_pattern.py -i $COMPLETEPATTERNFILE -k $ANNOTATIONFILE -p $PATTERN -o $OUTPUTFILE"
#echo "python $EXEDIR/find_gene_set_with_pattern_backbone.py -i $COMPLETEPATTERNFILE -k $ANNOTATIONFILE -p $PATTERN -n $N -o $OUTPUTFILE"