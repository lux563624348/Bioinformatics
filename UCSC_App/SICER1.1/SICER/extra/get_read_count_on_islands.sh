#!/bin/bash
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 7 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["Species"] ["FragmentSize"] ["ReadsFile"] ["IslandsFile"] ["OutputFile"] ["InputDir"] ["OutputDir"]
    echo ""
    exit 1
fi

SPECIES=$1
FRAGMENTSIZE=$2
READSFILE=$3
ISLANDSFILE=$4
OUTFILE=$5

INDIR=$6
OUTDIR=$7

echo "python $SICER/lib/associate_tags_with_regions.py -s $SPECIES -a $INDIR/$READSFILE -f $FRAGMENTSIZE  -b $INDIR/$ISLANDSFILE -o $OUTDIR/$OUTFILE"
python $SICER/lib/associate_tags_with_regions.py -s $SPECIES -a $INDIR/$READSFILE -f $FRAGMENTSIZE  -b $INDIR/$ISLANDSFILE -o $OUTDIR/$OUTFILE