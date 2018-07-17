#!/usr/bin/env python

"""
This code checks the correlation between the fc and the other features to see if there is any correlation. Limited to 2 LTR RLTR4s

assembled REs: 
	{chrom: 
		{(region_start, region_end):
			{"elements":[ids]; 
			"age":  
			target_name + "_rc":
			control_name + "_rc":
			"strand": "+";
			"num_boundary_elements": value; 
			"5_boundary_elements":[id]; 
			"5_boundary_elements_age":
			"3_boundary_elements":[id];
			"3_boundary_elements_age":
			"I_boundary_elements":[id]; 
			"expression_fc_" + target_name + "_vs_" + control_name: max_value
			}
		}
	}
	
86_derepressed_instances.pkl:{reClass:{reFamily:reName:[id]}}}

RE_instance_counts.pkl:{reClass:{reFamily:reName:total number of instances}}}
	

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import copy
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib
try:
   import cPickle as pickle
except:
   import pickle
import numpy as np
   
sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements')

import Utility_extended
import AnalyzeRNASeq
import age_comboRE_vs_derepression

def main(argv):
	#parser = OptionParser()
	#parser.add_option("-n", "--name_for_the_pickle_file", action="store", type="string", dest="summary_name", help="the name of the pickle file", metavar="<str>")
	#parser.add_option("-o", "--name_for_the_output_summary_file", action="store", type="string", dest="output_name", help="the name of the output file", metavar="<str>")
	#parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="RE tree in pickle format")
	
	#(opt, args) = parser.parse_args(argv)
	#if len(argv) < 6:
		#parser.print_help()
		#sys.exit(1)
		
	
	current_dir = os.getcwd()
	path = "/home/data/mm9/Lin/processed/RepElements/brokendown/summary"
	os.chdir(path)
	print "Loading the Assembled_REs.pkl"
	assembled_REs = pickle.load(open("Assembled_REs.pkl", 'rb'))
	num_re = 0
	for chrom in assembled_REs.keys():
		num_re += len(assembled_REs[chrom])
	print "There are %d elements in Assembled_REs.pkl" %num_re
	os.chdir(current_dir)
	
	num_LTR = 2
	two_LTR_assembled_REs = AnalyzeRNASeq.get_assembled_REs_subset_by_num_LTR(assembled_REs, num_LTR)
	
	
	feature_name = "5_boundary_elements_age"
	pc = 1
	title = feature_name + "_vs_" + "fc_on_assembled_RLTR4s" 
	list_feature, list_fc = AnalyzeRNASeq.correlate_feature_with_fc(two_LTR_assembled_REs, feature_name, pc)
	print len(list_fc)
	age_comboRE_vs_derepression.scatterplot(list_feature, list_fc, title, xscale='log', yscale='log')
	
	feature_name = "age"
	pc = 1
	title = feature_name + "_vs_" + "fc_on_assembled_RLTR4s" 
	list_feature_2, list_fc_2 = AnalyzeRNASeq.correlate_feature_with_fc(two_LTR_assembled_REs, feature_name, pc)
	print len(list_fc_2)
	age_comboRE_vs_derepression.scatterplot(list_feature_2, list_fc_2, title, xscale='log', yscale='log')
	
	feature_name = "5_boundary_elements_age"
	target_name = "index8"
	feature_name_b = target_name + "_rc"
	title = feature_name + "_vs_" + feature_name_b   + "_on_assembled_RLTR4s" 
	list_feature, list_feature_b = AnalyzeRNASeq.correlate_two_features(two_LTR_assembled_REs, feature_name, feature_name_b)
	age_comboRE_vs_derepression.scatterplot(list_feature, list_feature_b, title, xscale='log', yscale='log')
	
	feature_name = "age"
	feature_name_b = target_name + "_rc"
	title = feature_name + "_vs_" + feature_name_b   + "_on_assembled_RLTR4s" 
	list_feature, list_feature_b = AnalyzeRNASeq.correlate_two_features(two_LTR_assembled_REs, feature_name, feature_name_b)
	age_comboRE_vs_derepression.scatterplot(list_feature, list_feature_b, title, xscale='log', yscale='log')
	
	
if __name__ == "__main__":
	main(sys.argv)
