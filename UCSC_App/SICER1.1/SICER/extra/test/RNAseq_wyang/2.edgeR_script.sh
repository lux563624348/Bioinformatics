#! /bin/sh

file1=Set1_tags_on_Exons.txt

R --save<< EOT

library(edgeR)

raw.data <- read.delim("$file1.unique", header = FALSE)
matrix<-cbind(raw.data[,2],raw.data[,3])
f<-calcNormFactors(matrix, refColumn=1)
f<-f/exp(mean(log(f)))

D<-cbind(raw.data[,2],raw.data[,3])
group <- c("WT","KD2")
rownames(D) <- raw.data[, 1]
d<-DGEList(counts = D, group = group, lib.size = colSums(D)*f)
d<-estimateCommonDisp(d)
de.com <- exactTest(d)
write.table(topTags(de.com, n = 50000)$table, "Set1-de.txt", sep="\t")

EOT



