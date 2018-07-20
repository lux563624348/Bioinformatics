#!/usr/bin/env python
"""
This is a driver module to compare two lists of numbers to compare their distributions using MWU test
"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import scipy.stats

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
import gene_set_manipulation

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--datafile1", action="store", type="string", dest="datafile1", metavar="<file>", help="gene file 1, to be subtracted from")
	parser.add_option("-b", "--datafile2", action="store", type="string", dest="datafile2", metavar="<file>", help="gene file 2, to be subtracted")
	parser.add_option("-c", "--column1", action="store", type="int", dest="c1", metavar="<int>", help="the cloumn index for data file 1, 0-based")
	parser.add_option("-d", "--column2", action="store", type="int", dest="c2", metavar="<int>", help="the cloumn index for data file 2, 0-based")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	list1 = gene_set_manipulation.get_float_list(opt.datafile1, opt.c1)
	list2 = gene_set_manipulation.get_float_list(opt.datafile2, opt.c2)
	print "		p-value associated with ranking of the two sets of observations using MWU test: ", scipy.stats.mannwhitneyu(list1, list2)[1]
	
if __name__ == "__main__":
	main(sys.argv)