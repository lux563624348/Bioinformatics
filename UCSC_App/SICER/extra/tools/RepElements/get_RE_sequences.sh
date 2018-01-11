#!/bin/bash
EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RepElements
REDIR=/home/data/mm9/Lin/processed/RepElements/brokendown
GENOMEDIR=/home/data/mm9

RESUMMARYPICKLE=summary_on_LTR_ERV1_RLTR4_MM-int.pkl
GENOMEFILE=mm9.fa
OFILENAME=LTR_ERV1_RLTR4_MM-int.fa

python $EXEDIR/get_RE_sequences.py -i $REDIR/$RESUMMARYPICKLE -g $GENOMEDIR/$GENOMEFILE -o OFILENAME