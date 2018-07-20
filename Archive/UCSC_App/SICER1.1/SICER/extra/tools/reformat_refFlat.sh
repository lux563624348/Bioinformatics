#!/bin/bash

# reformat the refFlat format downloaded from the UCSC genome browser to the standard of our UCSC format.

#remove the comment line
sed 1d rmsk.txt > rmsk_comment_removed.txt

#move the gene symbol to the end
awk 'BEGIN { OFS = "\t"} { print $2,  $3, $4, $5, $6, $7, $8, $9, $10, $11, $1}' rmsk_comment_removed.txt > refFlat_mm10.ucsc

#sort the list according to the NM name
sort -g -k 1 refFlat_mm10.ucsc > sorted
mv sorted refFlat_mm10.ucsc

#Get the correspondance table between RefSeq name and gene symbol
cut -f1,11 refFlat_mm10.ucsc > refseq_vs_genesymbol_mm10.txt






