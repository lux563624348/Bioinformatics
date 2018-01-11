#!/bin/sh

#python /home/zang/Modules/normalize_to_1M_by_total.py -i $1-W$WINDOW_SIZE.graph -a 4 -t -1 -o normalized-$1-W$WINDOW_SIZE.graph 
sh variableStep.sh $1 $1.wig $2

####vstep for CHIP data:
#python /home/zang/Modules/normalize_to_1M_by_total.py -i SICER/$1.graph -a 4 -t -1 -o normalized-$1.graph
#sh variableStep.sh normalized-$1.graph normalized-$1.vstep $2

#sh variableStep.sh ../$1.graph $1.vstep $2
