#!/usr/bin/env python
"""
This is a driver module to overlap two sets of genes. Weiqun Peng
"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
import gene_set_manipulation

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--genefile1", action="store", type="string", dest="genefile1", metavar="<file>", help="gene file 1, to be subtracted from")
	parser.add_option("-b", "--genefile2", action="store", type="string", dest="genefile2", metavar="<file>", help="gene file 2, to be subtracted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="outputfile")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	list1 = gene_set_manipulation.get_gene_list(opt.genefile1, 0)
	list2 = gene_set_manipulation.get_gene_list(opt.genefile2, 0)
	result = gene_set_manipulation.gene_comparison(list1, list2)
	print opt.genefile1," number of unique ids: ", len(set(list1)) 
	print opt.genefile2," number of unique ids: ", len(set(list2))
	
	shared = result['shared']
	unique_in_1 = result['only in 1']
	unique_in_2 = result['only in 2']
	
	print "Number of shared ids are %d " %len(shared)
	print "Number of ids unique to %s are %d " %(opt.genefile1, len(unique_in_1))
	print "Number of ids unique to %s are %d " %(opt.genefile2, len(unique_in_2))
	
	#gene_set_manipulation.output_subset_in_file (opt.genefile1, opt.c1, IDs_to_be_retained, opt.outfile)
	
	
if __name__ == "__main__":
	main(sys.argv)