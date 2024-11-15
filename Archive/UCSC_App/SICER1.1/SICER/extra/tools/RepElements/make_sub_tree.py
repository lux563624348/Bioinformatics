#!/usr/bin/env python

"""

RE-tree is built by Characterize_RepElements.py

This module makes a subtree. The module takes two files are input:
The full RE_tree in pkl format
The file that is of the following tab delineated format:
Class	Family	Name

if Name is not there then  the whole family is in
if Family and Name are not there, then the who class is in.


"""


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import copy

try:
   import cPickle as pickle
except:
   import pickle

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')


sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements')

import UCSC_revised
import GenomeData
import SeparateByChrom
import Utility_extended
import RepElements
import associate_tags_with_regions
import get_total_tag_counts
import get_read_count_on_REs

plus = re.compile("\+");
minus = re.compile("\-");


def main(argv):
	parser = OptionParser()
	
	parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="file with RE tree in pickle format")
	parser.add_option("-l", "--RERetained", action="store", type="string", dest="REs_retained", metavar="<file>", help="file that contains the REs retained.")
	parser.add_option("-o", "--outfile", action="store", type="string",
			dest="outfile", help="outfile name", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	#load the RE tree to get the RE file names
	"""
	re_tree: {repClass:{repFamily:[repName]}}
	"""
	re_tree = pickle.load(open(opt.RE_Tree, 'rb'))
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	
	sub_tree = {}
	infile = open(opt.REs_retained, "r")
	for line in infile:
		line = line.strip();   #"""Strip white space off line""" 
		sline = line.split();
		if len(sline) == 1: # the whole class retained
			re_class = sline[0]
			assert (re_class in re_tree.keys())
			sub_tree[re_class] = copy.deepcopy(re_tree[re_class])
		elif len(sline) == 2: # the whole family retained
			re_class = sline[0]
			re_family = sline[1]
			assert (re_class in re_tree.keys())
			if re_class not in sub_tree.keys():
				sub_tree[re_class] ={}
			assert (re_family in re_tree[re_class].keys())
			sub_tree[re_class][re_family] = copy.deepcopy(re_tree[re_class][re_family])
		elif len(sline) == 3: # the whole family retained
			re_class = sline[0]
			re_family = sline[1]
			re_name = sline[2]
			assert (re_class in re_tree.keys())
			if re_class not in sub_tree.keys():
				sub_tree[re_class] ={}
			assert (re_family in re_tree[re_class].keys())
			if re_family not in sub_tree[re_class].keys():
				sub_tree[re_class][re_family] =[]
			sub_tree[re_class][re_family].append(re_name)
	
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(sub_tree)
	print "There are %d classes, %d family, and %d names in %s" %(numb_classes, numb_families, numb_names, opt.outfile)
	
	output = open(opt.outfile, 'wb')
	pickle.dump(sub_tree, output)
	output.close()	
	
if __name__ == "__main__":
	main(sys.argv)
