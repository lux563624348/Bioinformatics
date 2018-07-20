#!/bin/bash

EXEDIR=/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements

TREEDIR=/home/wpeng/data/mm9/Annotation/RepElements
TREE=mm9_re_tree.pkl

python $EXEDIR/make_sub_tree.py -t $TREEDIR/$TREE -l mysubtree.txt -o RLTR4_subtree.pkl
echo "python $EXEDIR/make_sub_tree.py -t $TREEDIR/$TREE -l mysubtree.txt -o RLTR4_subtree.pkl"