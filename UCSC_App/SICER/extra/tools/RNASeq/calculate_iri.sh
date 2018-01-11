#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq
EXE=calculate_iri.py
SUMMARYDIR=/home/data/hg19/JunZhu/processed/epigenome_roadmap/processed
SUMMARY=GSM669617_UCSF-UBC.CD4_Naive_Primary_Cells.mRNA-Seq.TC014_reformulated-on-EligibleEntrezGenes.dat_summary.pkl
NAME=CD4_Naive

echo "python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME-s hg19 "
python $EXEDIR/$EXE -f $SUMMARYDIR/$SUMMARY  -n $NAME -s hg19

