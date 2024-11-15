#!/bin/bash
RAWREADS=index6.sorted.bed


EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RepElements
TREEDIR=/home/data/mm9/Annotation/RepElements
REDIR=$TREEDIR/re_by_name
BEDDIR=/home/data/mm9/Lin/raw


SPECIES=mm9
UPSTREAMEXTENSION=0
DOWNSTREAMEXTENSION=0
FRAGMENTSIZE=0
TREE=mm9_re_tree.pkl


python $EXEDIR/get_read_count_on_RES.py -b $BEDDIR/$RAWREADS -f $FRAGMENTSIZE -t $TREEDIR/$TREE -l $REDIR -u $UPSTREAMEXTENSION -d $DOWNSTREAMEXTENSION -s $SPECIES 

	