#HiCcompare
library(HiCcompare)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

dat1 <- read.table(args[1], header=FALSE, col.names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "IF"))
dat2 <- read.table(args[2], header=FALSE, col.names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "IF"))

dat1 <- dat1[dat1$chr1==dat1$chr2 & dat1$chr1!="chrM" ]
dat2 <- dat2[dat2$chr1==dat2$chr2 & dat2$chr1!="chrM" ]

dat1 <- split(dat1, dat1$chr1)
dat2 <- split(dat2, dat2$chr1)

hic.list <- mapply(create.hic.table, dat1, dat2, SIMPLIFY = FALSE, scale=FALSE)

hic.list <- total_sum(hic.list)

register(MulticoreParam(workers = 10), default = TRUE)

hic.list <- hic_loess(hic.list, Plot=TRUE, parallel=TRUE)
hic.list <- hic_compare(hic.list, A.min = NA, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE, parallel=TRUE)

hic.list <- do.call(rbind, hic.list)

hic.list <- hic.list[hic.list$p.adj<0.05,]

write.table(hic.list, args[3], sep="\t", row.names=FALSE)
