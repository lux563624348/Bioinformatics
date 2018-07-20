#!/usr/bin/env python
"""
Assess the impact of mappability by collecting reads from all libraries.

Summary data structure
for each reClass_reFamily_reName:
{id:{feature_name:value}}

Given a summary name, print out its content in plain txt

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import copy

import numpy
import matplotlib.pyplot as plt
import matplotlib



try:
   import cPickle as pickle
except:
   import pickle


sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')

import AssembleFeatures
import AnalyzeRNASeq
import Utility_extended
		
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--name1", action="store", type="string", dest="summary_pkl_name", help="name for summary pickle file", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	assert( Utility_extended.fileExists(opt.summary_pkl_name) == 1)
	#features = ["index6_rc", "index8_rc"]
	#features = ["index12_rc", "index6_rc", "index8_rc", "index9_rc"]
	features = ["index12_rc", "index6_rc", "index8_rc", "index9_rc", "index3_rc"]
	
	
	summary = pickle.load(open(opt.summary_pkl_name, 'rb'))
	mykeys = sorted(summary.keys())
	total_rc = [0] * len(mykeys)
	for feature_name in features:
		feature_values = [ summary[myid][feature_name] for myid in mykeys]
		total_rc = [total_rc[i] + feature_values[i] for i in xrange(len(total_rc))]	
	log_total_rc = [log(i + 0.01, 10) for i in total_rc]
	plt.hist(log_total_rc, bins=100, color='r', alpha=1,label="Total Read Count" )	
	plt.savefig("Total_readcount_hist.png", format="png")
	plt.close()
	
	total_non_zero = 0
	for item in total_rc:
		if item >=1:
			total_non_zero += 1
	print "Total number of copies for this element is %d, total number of non-zero element is %d" %(len(total_rc), total_non_zero)
	
	
	summary_non_zero = {} #RE copies with non-zero read count in pooled data
	threshold = 5
	for i in xrange(len(total_rc)):
		if total_rc[i] >= threshold:
			summary_non_zero[mykeys[i]] = summary[mykeys[i]]
	
	AnalyzeRNASeq.get_log_fc_histogram(summary_non_zero, "index8", "index6", "RLTR4-mm-int-non-zero-rc.png", pc = 1)
	
	outfile_name = opt.summary_pkl_name + ".total_RC.dat"
	out = open(outfile_name, 'w')
	for i in xrange(len(total_rc)):
		if total_rc[i] >= 1:
			outw = str(mykeys[i]) + "\t" + str(total_rc[i]) + "\n"
			out.write(outw)
	out.close()
	
	
		
		
if __name__ == "__main__":
	main(sys.argv)
