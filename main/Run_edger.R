library(edgeR)
library("optparse")
#PATH='save/to/path'
#count_matrix = coverage # Colnames of coverage are sample names; each row of coverage represent a gene or peak
#groups = condition # A vector for each sample's corresponding condition

option_list <- list(
	make_option(c("-m","--input_matrix"),
		help="Please input a read count matrix"),
	make_option(c("-i","--group_index"),
		help="Please input a group label of matrix"),
	make_option(c("-o","--outpath"),
		help="Output File Path")
	)
opt <- parse_args(OptionParser(option_list=option_list))

print("#########################################################################")
print("Input Count Matrix File")
count_matrix <- opt$input_matrix
print (count_matrix)

print("Group File")
groups <- opt$group_index
print (groups)
print("#########################################################################")

print("Output File Path")
PATH <- opt$outpath
print (PATH)
print("#########################################################################")

experiment <- data.frame(exp=groups,row.names=colnames(count_matrix))
diff <- numeric(); N = length(unique(experiment$exp))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(experiment$exp==unique(experiment$exp)[i])
    condition2 = which(experiment$exp==unique(experiment$exp)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = count_matrix[,c(condition1,condition2)],group = group)
    cat(colnames(count_matrix[,condition1]));cat("  VS  ")
    cat(colnames(count_matrix[,condition2]));cat("\n")
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2.0))
    diff <- c(diff,diff_idx)
    change <- cbind(idx=diff_idx,
                    pval=round(qlf$table$PValue[diff_idx],7),
                    fdr=round(fdr[diff_idx],7),
                    logFC=round(qlf$table$logFC[diff_idx],7))
    write.table(change,paste(PATH,'/',unique(experiment$exp)[i],'-',unique(experiment$exp)[j],'.bed',sep=''),
               quote = F,sep = '\t',row.names = F,col.names = T)
  }
}
unique(diff)
