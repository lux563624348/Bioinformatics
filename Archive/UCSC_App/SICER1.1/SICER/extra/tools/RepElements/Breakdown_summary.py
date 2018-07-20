#!/usr/bin/env python

"""
Calculate read count on repetitive elements
Currently the data structure is {reClass:{reFamily:{reName:{id:{feature_name:value}}}}}, the memory usage is high and the speed is quite slow. With two features, the hard drive storage space is three times the storage used by {reClass:{reFamily:{reName:{id:value}}}}. The expected should be two times. The memory problem gets even worse when features are combined, 8G is not enough for combining index 6 and index 8. We need to do something. 
Solution: store read count for each type RE separately.
"""


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time

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

import Utility_extended
import get_read_count_on_REs


def main(argv):
	parser = OptionParser()
	
	parser.add_option("-i", "--summary_pickle_file", action="store", type="string", dest="summary", metavar="<file>", help="summary file in pickle format")
	parser.add_option("-n", "--lib_name", action="store", type="string", dest="name",help="name of the library", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	if Utility_extended.fileExists(opt.summary):
		my_summary = pickle.load(open(opt.summary, 'rb'))
		(numb_classes, numb_families, numb_names, numb_ids, features) = get_read_count_on_REs.characteristics(my_summary)
		print "There are %d re classes, %d re family, and %d re names in %s." %(numb_classes, numb_families, numb_names, opt.summary)
		print "There are %d re instances in %s." %(numb_ids, opt.summary)
		print "There are %d features in summary: %s" %(len(features), ",".join(features))
	
		get_read_count_on_REs.breakdown_and_output(my_summary, opt.name)
		
	else:
		print opt.summary, " is not found";
		sys.exit(1)	
	
if __name__ == "__main__":
	main(sys.argv)
