#!/bin/bash
# 
PATHTO=/home/data/SICER1.1
SICER=$PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

LIB1=$1
LIB2=$2
COL1=$3
COL2=$4

echo "python $SICER/extra/tools/compare_two_lists.py -a $LIB1 -b $LIB2 -c $COL1 -d $COL2"
python $SICER/extra/tools/compare_two_lists.py -a $LIB1 -b $LIB2 -c $COL1 -d $COL2