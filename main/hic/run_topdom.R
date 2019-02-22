#!Given some moron's contribution, this file is necessary.
library("optparse")

option_list <- list(
	make_option(c("-i","--input_matrix"),
		help="Please input a nx(n+3) HiC interaction matrix"),
	make_option(c("-w","--window_size"),default=5,metavar="number",type="integer",
		help="Please choose window size between 5~20"),
	make_option(c("-o","--outfile"),
		help="Output File Path")
	)
opt <- parse_args(OptionParser(option_list=option_list))



print("#########################################################################")
print("Input Matrix File")
print (opt$input_matrix)
print("Windows Size")
print (opt$window_size)
print("#########################################################################")


source("TopDom_v0.0.2.R")
TopDom(matrix.file = opt$input_matrix, window.size= opt$window_size, outFile=opt$outfile)

