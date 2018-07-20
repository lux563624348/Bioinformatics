#!/bin/bash
# 
# Authors: Sean Grullon,  Weiqun Peng
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


# The path to the data files.
DATADIR=$1

# Input sample bed file
RAWBED=$2

# Species, for allowed species see GenomeData.py
SPECIES=$3

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=$4

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
echo " "
echo " "
echo "Automatically calculate window size ..."
echo "python $SICER/utility/calculate-window-size.py -s $SPECIES -b $DATADIR/$RAWBED -f $FRAGMENT_SIZE"
let WINDOW_SIZE=`python $SICER/utility/calculate-window-size.py -s $SPECIES -b $DATADIR/$RAWBED -f $FRAGMENT_SIZE | grep Window | awk '{print $2}'`
#let GAP_SIZE=$WINDOW_SIZE*$GAP_SIZE

echo " " 
echo "Our Window Size is: $WINDOW_SIZE"
