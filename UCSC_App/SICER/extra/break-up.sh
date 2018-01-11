#!/bin/sh
for ((i=1;i<=22;i+=1)); do
    grep chr$i[[:space:]] ../CD4-IgG.bed > chr$i-IgG.bed
done

grep chrX[[:space:]] ../CD4-IgG.bed > chrX-IgG.bed
grep chrY[[:space:]] ../CD4-IgG.bed > chrY-IgG.bed

