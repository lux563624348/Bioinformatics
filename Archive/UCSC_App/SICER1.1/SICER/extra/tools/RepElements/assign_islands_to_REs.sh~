#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RepElements

UPSTREAM=100
DOWNSTREAM=100
SPECIES=mm9
TREEDIR=/home/data/hg19/Annotation/RepElements
TREE=mm9_re_tree.pkl
REDIR=$TREEDIR/re_by_name

echo "python $EXEDIR/assign_islands_to_REs -b $BEDDIR/$ISLAND -t $TREEDIR/$TREE -l $REDIR -u $UPSTREAM -d $DOWNSTREAM -s $SPECIES"
python $EXEDIR/assign_islands_to_REs -b $BEDDIR/$ISLAND -t $TREEDIR/$TREE -l $REDIR -u $UPSTREAM -d $DOWNSTREAM -s $SPECIES