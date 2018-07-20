#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/pattern_analysis
INFILESUFFIX=-Promoter-ReadCount.dat
FEATURELIST=featurelist.txt

python $EXEDIR/make-01-pattern-data.py -i $INFILESUFFIX -l $FEATURELIST -o PromoterDigitalPattern.txt

echo "python $EXEDIR/make-01-pattern-data.py -i $INFILESUFFIX -l $FEATURELIST -o PromoterDigitalPattern.txt"