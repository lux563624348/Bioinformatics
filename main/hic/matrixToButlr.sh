#!/bin/bash

EXEDIR=/home/zzeng/Software/BUTLRTools-master

echo "perl $EXEDIR/matrixToButlr.pl -g mm9_genome_size.txt -m WT_CD8_intrachrom_matrix_list.txt -a mm9 -r 10kb -o WT_CD8.10kb.btr"
perl $EXEDIR/matrixToButlr.pl -g mm9_genome_size.txt -m WT_CD8_intrachrom_matrix_list.txt -a mm9 -r 10kb -o WT_CD8.10kb.btr
echo ""
echo ""

#echo "perl $EXEDIR/matrixToButlr.pl -g mm9_genome_size.txt -m DKO_CD8_intrachrom_matrix_list.txt -a mm9 -r 10kb -o DKO_CD8.10kb.btr"
#perl $EXEDIR/matrixToButlr.pl -g mm9_genome_size.txt -m DKO_CD8_intrachrom_matrix_list.txt -a mm9 -r 10kb -o DKO_CD8.10kb.btr
#echo ""
#echo ""

echo "done"
