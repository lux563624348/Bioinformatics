#!/usr/bin/env python

import pandas as pd
from optparse import OptionParser

def BED_Correction(INPUT_BED, GENOME_SIZE_FILE):
    df = pd.read_csv(INPUT_BED, sep='\t',header=None)
    df_size = pd.read_csv(GENOME_SIZE_FILE, sep='\t',header=None)
    for chromosome in df[0].unique():
        genome_length = df_size[df_size[0]==chromosome][1].values
        max_length = df.iloc[df[df[0]==chromosome].tail(1).index, 2].values
        if (max_length > genome_length):
            df.iloc[df[df[0]==chromosome].tail(1).index, 2] = genome_length
    df.to_csv(INPUT_BED[:-4:]+ ".track.bed", sep='\t', header=None, index=None)
    return "BED Corretion For " + INPUT_BED[:-4:]


def main(argv):
	desc="Bed Correction"
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="input_path", metavar="<file>")
	parser.add_option("-s", "--genome_size", action="store", type="string",
		dest="genome_size", help="genome_size", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	
	print " "
	print "Here is the Summary of your input."
	print "Input Path Contains BED need to be changed: %s" % opt.input_path
	print "Genomesize file: %s" % opt.genome_size
	
	BED_Correction(opt.input_path, opt.genome_size)




if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
