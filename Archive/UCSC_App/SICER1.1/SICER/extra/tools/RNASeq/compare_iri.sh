#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq
EXE=compare_iri.py
SUMMARYDIR=/home/guestuser/IRI_analysis/lip_Leonard
SUMMARY=hCD8_4h.txt
NAME=hCD8_4h

echo "python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME -s hg19 "
python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME -s hg19


SUMMARY=hCD8_24h.txt
NAME=hCD8_24h

echo "python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME -s hg19 "
python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME -s hg19
