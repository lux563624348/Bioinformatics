#!/bin/bash
# Authors: Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).


##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

EXEDIR=/home/data/SICER1.1/SICER/extra
EXE=fft.py
OUTFILE=$1-rfft.dat

echo "python $EXEDIR/$EXE $1 $OUTFILE"

python $EXEDIR/$EXE $1 $OUTFILE