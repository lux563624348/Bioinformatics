#!/bin/bash
SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib:$SICER/extra
export PYTHONPATH
echo ""
echo ""
echo "Bivalent domain Promoters in WT"
PATTERN=11nn
sh find_geneset_w_pattern.sh $PATTERN


echo "Bivalent domain lost K4 in KO"
PATTERN=110n
sh find_geneset_w_pattern.sh $PATTERN

echo "Bivalent domain lost K9 in KO"
PATTERN=11n0
sh find_geneset_w_pattern.sh $PATTERN


