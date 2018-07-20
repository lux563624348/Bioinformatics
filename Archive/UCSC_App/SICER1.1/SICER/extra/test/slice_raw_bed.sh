#!/bin/bash
# Copyright (c) 2011 The George Washington University
# Authors: Weiqun Peng
#

SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 3 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["Raw read file "] ["Number of tags needed"] {"Output read file"]
    echo ""
    exit 1
fi

SAMPLEBED=$1
SAMPLE=${SAMPLEBED%.*}

echo "python $SICER/utility/slice_raw_bed.py -f $SAMPLEBED -n $2 -o $3"
python $SICER/utility/slice_raw_bed.py -f $SAMPLEBED -n $2 -o $3