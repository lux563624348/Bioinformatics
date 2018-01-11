#!/usr/bin/env python
"""
data structure
for each reClass_reFamily_reName:
{id:{feature_name:value}}

Assemblage of annotation into the summary is done in AnalyzeRNASeq.py

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import copy

try:
   import cPickle as pickle
except:
   import pickle


sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')


import Utility_extended
import get_read_count_on_REs

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")

def combine_features(new_feature, summary = None):
	"""
	new_feature: {id:{feature_name:value}}
	summary: {id:{feature_name:value}}
	"""
	if summary is None:
		return new_feature

	else:
		for myid in new_feature.keys():
			summary[myid].update (new_feature[myid])
		return summary

def print_out_summary (lib_file, outputfile):
	"""
	lib file in pickle format: {repID:{feature:value}}
	output the rc data of the particular pickle file 
	"""
	assert( Utility_extended.fileExists(lib_file) == 1)
	lib = pickle.load(open(lib_file, 'rb'))
	myid = lib.keys()[0]
	mykeys = (lib[myid]).keys()
	mykeys = sorted(mykeys)
	of = open (outputfile, 'w')
	oline = "ID" + "\t" + ("\t").join(mykeys) + "\n"
	of.write(oline)
	for myid in lib.keys():
		re  = lib[myid]
		oline = str(myid)	
		for feature in mykeys:
			oline +=  "\t" + str(re[feature])
		oline += "\n"
		of.write(oline)
	of.close()		
		
		
def print_out(name, reClass, reFamily, reName, outputfile):
	"""
	{repID:{feature:value}}
	output the rc data of particular (repClass, repFamily, repName)
	"""
	lib_file = name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
	print_out_summary (lib_file, outputfile)		
		
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--name1", action="store", type="string", dest="summary_pkl_name", help="name for combined result pickle file", metavar="<file>")
	parser.add_option("-j", "--name2", action="store", type="string", dest="new_lib_name", help="name for additional summary pickle file", metavar="<file>")
	parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="RE tree in pickle format")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	#load the RE tree to get the RE file names
	re_tree = pickle.load(open(opt.RE_Tree, 'rb'))
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	
	
	for reClass in re_tree.keys():
		for reFamily in re_tree[reClass].keys():
			for reName in re_tree[reClass][reFamily]:
				summary_file_name = opt.summary_pkl_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				if Utility_extended.fileExists(summary_file_name) != 1:
					print summary_file_name
					exit(1)
				summary = pickle.load(open(summary_file_name, 'rb'))
				
				new_lib_file_name = opt.new_lib_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(new_lib_file_name) == 1)
				new_lib = pickle.load(open(new_lib_file_name, 'rb'))
				
				#new_feature: {id:{feature_name:value}}
				#summary: {id:{feature_name:value}}
				summary = combine_features(new_lib, summary)
				outfile_name = "summary" + "_on_"  + "_".join([reClass, reFamily, reName]) + ".pkl"
				output = open(outfile_name, 'wb')
				pickle.dump(summary, output)
				output.close()	

	repClass =re_tree.keys()[0]
	repFamily = re_tree[repClass].keys()[0]
	repName = re_tree[repClass][repFamily][0]
	outfile_name = "summary_on_" +  "_".join([repClass, repFamily, repName]) + ".dat"
	print_out("summary", repClass, repFamily, repName, outfile_name)

if __name__ == "__main__":
	main(sys.argv)
