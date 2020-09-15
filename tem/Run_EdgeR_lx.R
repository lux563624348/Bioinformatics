library("edgeR")
library("optparse")
## input file format is the same as edger from Galaxy.

option_list <- list(
  make_option(c("-c","--count_table"),
              help="Please input a read count matrix"),
  make_option(c("-i","--group_index"),
              help="Please input a group label of matrix"),
  make_option(c("-o","--outpath"),
              help="Output File Path")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("#########################################################################")
print("Input File of Count Table")
PATH_count_table=opt$count_table
print (PATH_count_table)

print("Input Index File")
PATH_group_index <- opt$group_index
print (PATH_group_index)
print("#########################################################################")

print("Output File Path")
Out_PATH <- opt$outpath
print (Out_PATH)
print("#########################################################################")


rawCountTable = read.table(PATH_count_table, header=TRUE, sep="\t", row.names=1)
sampleInfo = read.table(PATH_group_index, header=TRUE, sep="\t", row.names=1)

Out_PATH=Out_PATH


group = factor(sampleInfo$Genotype)
y = DGEList(counts=rawCountTable,group=group)
y_status = filterByExpr(y)
y = calcNormFactors(y)
design = model.matrix(~group)
y = estimateDisp(y,design)
fit = glmQLFit(y,design)
qlf = glmQLFTest(fit,coef=2)
FDR = p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))

output_table = cbind(GeneID=rownames(qlf$table), logFC=round(qlf$table$logFC,4),
        PValue=round(qlf$table$PValue,6),FDR=round(FDR,6),Counts_Pass=y_status)


write.table(output_table, paste(Out_PATH,'/', unique(sampleInfo$Genotype)[2],"_vs_",unique(sampleInfo$Genotype)[1],".txt", sep="" ), quote=F, sep='\t', row.names = F, col.names = T)