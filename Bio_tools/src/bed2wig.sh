#!/bin/bash

#1. generate graph file from raw bed:
SICER=~/Software/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

SAMPLEDIR=$1
SAMPLE=$2
SAMPLEBED=$SAMPLE.bed
WINDOW_SIZE=$3
FRAGMENT_SIZE=$4
SPECIES=$5
SUMMARY=$SAMPLE-W$WINDOW_SIZE.graph
SUMMARYWIG=$SAMPLE-W$WINDOW_SIZE.wig
NORMALIZEDSUMMARY=$SAMPLE-W$WINDOW_SIZE-normalized.graph
NORMALIZEDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-RPKM.wig

cd $SAMPLEDIR

echo "Partion the genome in windows and generate summary files ..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY "
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY 

# echo ""
# echo ""
# echo "Convert the sample summary graph into wig vstep format..."
# echo "sh $SICER/src/variableStep.sh $SUMMARY $SUMMARYWIG $SAMPLE $WINDOW_SIZE"
# sh $SICER/src/variableStep.sh $SUMMARY $SUMMARYWIG $SAMPLE $WINDOW_SIZE


scaling=$( echo "1000000000.0/$WINDOW_SIZE" | bc)


echo ""
echo ""
echo "Normalize summary graph to RPKM for $SAMPLEBED ..."
echo "python $SICER/src/normalize.py -i $SUMMARY -a 3 -t $scaling -o $NORMALIZEDSUMMARY"
python $SICER/src/normalize.py -i $SUMMARY -a 3 -t $scaling -o $NORMALIZEDSUMMARY



echo ""
echo ""
echo "Convert the normalized sample summary graph into wig vstep format..."
echo "sh $SICER/src/variableStep.sh $NORMALIZEDSUMMARY $NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $NORMALIZEDSUMMARY $NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE


rm $NORMALIZEDSUMMARY
rm $SUMMARY


