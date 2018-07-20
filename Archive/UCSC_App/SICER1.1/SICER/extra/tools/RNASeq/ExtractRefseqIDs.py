#!/usr/bin/env python



import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
from operator import itemgetter
import scipy.stats
import numpy
import matplotlib.pyplot as plt

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import AnalyzeRNASeqPASeq

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--entrez_genes", action="store", type="string", dest="one", help="combined result file", metavar="<file>")
	parser.add_option("-p", "--option", action="store", type="string", dest="option", help="main or all", metavar="<str>")
	parser.add_option("-o", "--Refseq IDs", action="store", type="string", dest="Refseq_ID_file", help="refseq ID file", metavar="<file>")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	One = AnalyzeRNASeqPASeq.readfile(opt.one)
	if opt.option == "main":
		AnalyzeRNASeqPASeq.getMainRefseqIDs(One, opt.Refseq_ID_file)
	elif opt.option == "all":
		AnalyzeRNASeqPASeq.getAllRefseqIDs(One, opt.Refseq_ID_file)
	else:
		print opt.option, " can only be main or all"
		exit(1)

if __name__ == "__main__":
	main(sys.argv)