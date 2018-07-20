#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/Repeats
DATADIR=/home/data/hg19/Annotation/RepElements

python $EXEDIR/RepElements.py -i $DATADIR/RepElements_hg19.txt -s hg19 -o temp