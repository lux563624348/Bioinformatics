#!/bin/bash

awk 'BEGIN { OFS = "\t"} { print $2,  $3, $4, $5, $6, $7, $8, $9, $10, $11, $1}' refFlat_hg19.txt > refFlat_hg19.ucsc

cut -f3,7 refLink_hg19.txt > refseq_vs_EntrezID_hg19.txt

sort -k 1 refFlat_hg19.ucsc > refFlat_hg19_sorted.ucsc
sort -k 1 refseq_vs_EntrezID_hg19.txt > refseq_vs_EntrezID_hg19_sorted.txt


# why is --nocheck-order is needed is not clear
# => Now I found out: it is because one refseq ID might have multiple entries in the file
# Once both files are sorted, join is able to ignore lines unique in one file.
# if an id appears multiple times in one file, and one times in the other file, the output will keep the multiple 
# appearance of the id.

join  --nocheck-order refFlat_hg19_sorted.ucsc refseq_vs_EntrezID_hg19_sorted.txt > refFlat_hg19_EntrezID.ucsc

exedir=/home/data/SICER1.1/SICER/extra

Due to concern about join, I've implemented filter_refFlat.py that
1) remove the non-canonical chrom entries
2) remove redundancy in refseq entries
3) join the refFlat and Entrez ID

The end result is refFlat_hg19_EntrezID_filtered.ucsc


