#!/usr/bin/env python
"""
This is a driver module to subtract a subset. Weiqun Peng
"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
import Utility_extended

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="infile", metavar="<file>", help="name for input file, which is a multicolumn text file")
	parser.add_option("-c", "--column", action="store", type="int", dest="c", metavar="<int>", help="the index of the column to be rescaled, 1-based")
	parser.add_option("-r", "--rescale_factor", action="store", type="float", dest="rescale_factor", metavar="<float>", help="the rescale factor that will be multiplied to the numbers in column c")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="name for output file")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	Utility_extended.rescale_a_column(opt.infile, opt.c-1, opt.rescale_factor, opt.outfile)
	
if __name__ == "__main__":
	main(sys.argv)