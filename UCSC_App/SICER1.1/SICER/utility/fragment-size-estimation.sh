#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

PATHTO=/SICER1.1
SICER=$PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 2 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["bed file"] ["species"]
    echo ""
    exit 1
fi

# INPUT is the sample bed file
INPUT=$1
SAMPLE=${INPUT%.*}

SPECIES=$2

# Directory where the INPUT bed file is. local as default
DIR=.

# RESOLUTION is the data resolution of the estimated fragment size (in bp). It is the window size of the summary graph file generated and the step size of the cross correlation calculation. Decreasing this number would significantly increase the computing time.
RESOLUTION=10


echo "########################################################################################"
echo "##  Estimating ChIP fragment size by cross correlation between watson and crick tags  ##"
echo "##                          An affiliated tool of SICER v1.1                          ##"
echo "########################################################################################"
echo "Data resolution: $RESOLUTION"

# Calculating tag correlation in the longest chromosome in the given species

CHROM=`python $SICER/utility/lookup_chrom_length.py -s $SPECIES | grep $SPECIES | awk '{print $2}'`
let CHROM_LENGTH=`python $SICER/utility/lookup_chrom_length.py -s $SPECIES | grep $SPECIES | awk '{print $3}'`
echo ""
echo "The longest chromosome in $SPECIES is $CHROM with length $CHROM_LENGTH"
echo ""

# Output file name
DATAFILE=$SAMPLE-$CHROM-tag-correlation

# Grep tags for longest chromosome
echo "Separate $CHROM tags to calculate fragment size ..."
grep ${CHROM}[[:space:]] $DIR/$INPUT > temp1.bed
echo ""
echo "Total tag count in $CHROM is"
wc -l temp1.bed
echo ""

# Separate watson and crick reads
echo "Separate watson and crick reads in $CHROM ..."
grep [[:space:]]+ temp1.bed > temp_plus.bed
grep [[:space:]]- temp1.bed > temp_minus.bed
echo ""

# Make summary graph files
echo "Make summary graph files"
echo ""
echo "Positive"
echo "python $SICER/lib/make_graph_file.py -f temp_plus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_plus.graph"
python $SICER/lib/make_graph_file.py -f temp_plus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_plus.graph
echo ""

echo ""
echo "Negative"
echo "python $SICER/lib/make_graph_file.py -f temp_minus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_minus.graph"
python $SICER/lib/make_graph_file.py -f temp_minus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_minus.graph
echo ""

# Calculate cross correlation between watson summary graph and crick summary graph
echo ""
echo "Calculating correlation between watson and crick summary graph for fragment size estimation"
echo "python $SICER/utility/calculate_cross_correlation_long_range.py -s $SPECIES -c $CHROM -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $DIR/$DATAFILE.txt"
let FRAGMENT_SIZE=`python $SICER/utility/calculate_cross_correlation_long_range.py -s $SPECIES -c $CHROM -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $DIR/$DATAFILE.txt | grep fragment | awk '{print $2}'`

echo "Done! The estimated ChIP fragment size is $FRAGMENT_SIZE. Please double check by reading the peak location in the plot ${DATAFILE}_plot.eps"


rm temp1.bed
rm temp_minus.bed
rm temp_minus.graph
rm temp_plus.bed
rm temp_plus.graph
