#!/bin/bash
# Copyright (c) 2011 The George Washington University
# Authors: Weiqun Peng
#

DIR=/home/data/hg18/CD16/raw/set1
n=500000

for INPUT in GA1564_hg18_CD16-H3K4me1-A GA1566_hg18_CD16-IgG-A GA1568_hg18_CD16-H3K4me1-B GA1570_hg18_CD16-IgG-B
do
sh  slice_raw_bed.sh $DIR/$INPUT.bed $n $INPUT-n$n.bed >> log
echo "sh  slice_raw_bed.sh $DIR/$INPUT.bed $n $INPUT-n$n.bed>> log"
echo ""
echo ""
echo ""
echo ""

done
