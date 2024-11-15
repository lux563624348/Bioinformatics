#!/usr/bin/env python

"""
ChIP-Seq can be multiple ChIP-Seq libraries:
The 


summary pickle data structure 
{id:
	{'annotation':RepElement class instance
	'H3K4me3_WT-W200-G200_rc':
	'H3K4me3_WT-W200-G200_rpkm':
	'H3K4me3_mir34bc_KO-W200-G200_rc':
	'H3K4me3_mir34bc_KO-W200-G200_rpkm':
	'H3K9me3_WT-W200-G400_rc':
	'H3K9me3_WT-W200-G400_rpkm':
	'H3K9me3_mir34bc_KO-W200-G400_rc':
	'H3K9me3_mir34bc_KO-W200-G400_rpkm':	
	}
}

output:

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
import RepElements
import AssembleFeatures
import get_read_count_on_REs
import AnalyzeRNASeq


def calculate_level(reClass_reFamily_reName_summary, feature_name):
	"""
	reClass_reFamily_reName_summary: {id:{feature_name:value}}
	return (average, median)
	"""
	mylist = [ reClass_reFamily_reName_summary[myid][feature_name] for myid in reClass_reFamily_reName_summary.keys()]
	return np.average(mylist), np.median(mylist)
	
def calculate_enrichment(reClass_reFamily_reName_summary, feature_name, min_level=0.00001):
	"""
	enrichment_ratio = # sites with the mark/ # of sites
	"""
	mylist = [ reClass_reFamily_reName_summary[myid][feature_name] for myid in reClass_reFamily_reName_summary.keys()]
	total = len(mylist)
	num_above_min = 0
	for item in mylist:
		if item > min_level:
			num_above_min += 1
	return num_above_min*1.0/total
	
def get_feature_level(re_tree, summary_name):
	"""
	Find the mean rpkm, median rpkm, or presence of each mark for each species 
	return:
	feature_level_mean: {reClass:{reFamily:{reName:{feature_name:level}}}}
	feature_level_median: {reClass:{reFamily:{reName:{feature_name:level}}}}
	feature_enrichment: {reClass:{reFamily:{reName:{feature_name:enrichment_ratio}}}}, 
		where enrichment_ratio = # sites with the mark/ # of sites
	"""
	feature_level_mean = {}
	feature_level_median = {}
	feature_enrichment = {}
	
	flag = 0
	
	for reClass in re_tree.keys():
		feature_level_mean[reClass] = {}
		feature_level_median[reClass] = {}
		feature_enrichment[reClass] = {}
		
		for reFamily in re_tree[reClass].keys():
			feature_level_mean[reClass][reFamily] = {}
			feature_level_median[reClass][reFamily] = {}
			feature_enrichment[reClass][reFamily] = {}
			
			for reName in re_tree[reClass][reFamily]:
				feature_level_mean[reClass][reFamily][reName] = {}
				feature_level_median[reClass][reFamily][reName] = {}
				feature_enrichment[reClass][reFamily][reName] = {}
				
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				# {id:{feature_name:value}}
				reClass_reFamily_reName_summary = pickle.load(inf)
				inf.close()
				
				if flag == 0: # Do this only one time
					
					feature_names = AnalyzeRNASeq.get_feature_names(reClass_reFamily_reName_summary)
					print "\nFeature names are: ",  feature_names
					flag = 1
					
				for feature_name in feature_names:
					if feature_name != "annotation":
						mean, median = calculate_level(reClass_reFamily_reName_summary, feature_name)
						feature_level_mean[reClass][reFamily][reName][feature_name] = mean
						feature_level_median[reClass][reFamily][reName][feature_name] = median
						feature_enrichment[reClass][reFamily][reName][feature_name] = calculate_enrichment(reClass_reFamily_reName_summary, feature_name)
					
				
	return feature_level_mean, feature_level_median, feature_enrichment
	
	
def rank_by_feature_level (feature_level_dic, feature_name, out_file_name):
	"""
	feature_level_dic: {reClass:{reFamily:{reName:{feature_name:level}}}}
	"""
	mylist = []
	for reClass in feature_level_dic.keys():
		for reFamily in feature_level_dic[reClass].keys():
			for reName in feature_level_dic[reClass][reFamily]:
				level = feature_level_dic[reClass][reFamily][reName][feature_name]
				mylist.append((reClass, reFamily, reName, level))
	mylist.sort(key = itemgetter(3), reverse=True)
	if out_file_name != "":
		outf = open(out_file_name, "w")
		for item in mylist:
			outf.write("\t".join(map(str,item)) + "\n")
		outf.close()
	return mylist
	
def find_pattern(re_tree, summary_name, present_list, absent_list, threshold = 0.0001):
	"""
	present_list:[feature_name] features that are required to be present
	absent_list:[feature_name] features that are required to be absent
	
	return: {reClass:{reFamily:{reName:{feature_name:enrichment_ratio}}}}
	{reClass:{reFamily:{reName:{feature_name:[ids]}}}}
	"""
	
	pattern_enrichment = {}
	pattern_positve_ids = {}
	flag = 0
	
	for reClass in re_tree.keys():
		pattern_enrichment[reClass] = {}
		pattern_positve_ids[reClass] = {}
		
		for reFamily in re_tree[reClass].keys():
			pattern_enrichment[reClass][reFamily] = {}
			pattern_positve_ids[reClass][reFamily] = {}
			
			for reName in re_tree[reClass][reFamily]:
				
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				# {id:{feature_name:value}}
				reClass_reFamily_reName_summary = pickle.load(inf)
				inf.close()
				
				if flag == 0: # Do this only one time
					
					feature_names = AnalyzeRNASeq.get_feature_names(reClass_reFamily_reName_summary)
					print "\nFeature names are: ",  feature_names
					
					assert (set(present_list).issubset(set(feature_names)))
					assert (set(absent_list).issubset(set(feature_names)))
					
					flag = 1
				
				enrichment, positive_ids = calculate_pattern_enrichment_in_single_species(reClass_reFamily_reName_summary, present_list, absent_list, threshold)
				
				pattern_enrichment[reClass][reFamily][reName] = enrichment
				pattern_positve_ids[reClass][reFamily][reName] = positive_ids
	
	return pattern_enrichment, pattern_positve_ids
	
						
def calculate_pattern_enrichment_in_single_species(reClass_reFamily_reName_summary, present_list, absent_list, threshold = 0.0001):
	"""
	Find the number of occurances in which marks in present_list are present and in absent_list are absent.
	# reClass_reFamily_reName_summary: {id:{feature_name:value}}
	"""
	negative_list = []
	total = len(reClass_reFamily_reName_summary.keys())
	total_negative = 0
	for myid in reClass_reFamily_reName_summary.keys():
		negative = 0
		for feature_name in present_list:
			if reClass_reFamily_reName_summary[myid][feature_name] <= threshold:
				negative = 1
				negative_list.append(myid)
				total_negative += 1
				break
		if negative == 0:
			for feature_name in absent_list:
				if reClass_reFamily_reName_summary[myid][feature_name] > threshold:
					negative_list.append(myid)
					total_negative += 1
					break
	assert (len(negative_list) == total_negative)
	
	occurance = total-total_negative
	positive_list = list(set(reClass_reFamily_reName_summary.keys()) - set(negative_list))
	enrichment = occurance*1.0/total
	return enrichment, positive_list
	
def rank_by_pattern_enrichment(pattern_enrichment, out_file_name):
	
	"""
	pattern_enrichment: {reClass:{reFamily:{reName:enrichment}}}
	mylist: [(reClass, reFamily, reName, enrichment)] sorted by enrichment
	"""
	mylist = []
	for reClass in pattern_enrichment.keys():
		for reFamily in pattern_enrichment[reClass].keys():
			for reName in pattern_enrichment[reClass][reFamily]:
				enrichment = pattern_enrichment[reClass][reFamily][reName]
				mylist.append((reClass, reFamily, reName, enrichment))
	mylist.sort(key = itemgetter(3), reverse=True)
	if out_file_name != "":
		outf = open(out_file_name, "w")
		for item in mylist:
			outf.write("\t".join(map(str,item)) + "\n")
		outf.close()
	return mylist

def output_pattern_ids(pattern_positve_ids, reClass, reFamily, reName, out_filename):
	"""
	pattern_positve_ids:{reClass:{reFamily:{reName:{feature_name:[ids]}}}}
	"""
	outf = open(out_filename, "w")
	for myid in pattern_positve_ids[reClass][reFamily][reName]:
		outline = myid + "\n"
		outf.write(outline)
	outf.close()

	
def filter_ranked_list_by_RE_instance_count(mylist, RE_instance_counts, out_file_name, min_count=0):
	"""
	mylist: [(reClass, reFamily, reName, enrichment)]
	RE_instance_counts: {reClass:{reFamily:{reName:counts}}}
	
	return filtered list
	"""
	filtered_list = []
	for item in mylist:
		reClass = item[0]
		reFamily = item[1]
		reName = item[2]
		enrichment = item[3]
		if RE_instance_counts[reClass][reFamily][reName] >= min_count:
			filtered_list.append((reClass, reFamily, reName, enrichment, RE_instance_counts[reClass][reFamily][reName]))
			
	outf = open(out_file_name, "w")
	for item in filtered_list:
		outf.write("\t".join(map(str,item)) + "\n")
	outf.close()
	return filtered_list
		
def plot_top_N(mylist, number, name):
	"""
	Plot bar plot: x: name, y:enrichment
	mylist: [(reClass, reFamily, reName, enrichment)]
	"""
	OY = [item[3] for item in mylist[:number]]
	OX = [("-").join(item[:3]) for item in mylist[:number]]
	fig = plt.figure()
	width = .35
	ind = np.arange(len(OY))
	plt.bar(ind, OY)
	plt.xticks(ind + width / 2, OX)
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	#fig.autofmt_xdate()
	plt.savefig(name)
	
	
						
def find_differentially_marked_families(target_name, control_name, suffix, re_tree, summary_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc, pc = 1):
	"""
	calculate the de-enrichment for each re family and then rank according to enrichment, output the top10
	
	"""
	
	print "\n\n", target_name, " vs ",  control_name
	full_target_name = target_name + suffix
	full_control_name = control_name + suffix
	
	(de_members, de_enrichment, total_instances) = AnalyzeRNASeq.get_all_upregulated_REs(re_tree, summary_name, full_target_name, full_control_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc, pc)
	
	out_file_name = target_name + "_vs_" + control_name + "_fc" + str(fc) + "DE_enrichment_ranking.txt"
	#mylist: [(reClass, reFamily, reName, de_enrichment[reClass][reFamily][reName], total_instances[reClass][reFamily][reName])]
	mylist = AnalyzeRNASeq.rank_enrichment(de_members, de_enrichment, total_instances, out_file_name)
	
	
	for i in xrange(min(10,len(mylist))):
		reClass = mylist[i][0]
		reFamily = mylist[i][1]
		reName = mylist[i][2]
		number_instances = mylist[i][4]
	
		if number_instances > 50: 
			summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
			inf = open(summary_file_name, 'rb')
			summary = pickle.load(inf)
			if i == 0:
				print AnalyzeRNASeq.get_feature_names(summary)
			inf.close()
			
			# output instances of differentially expressed regions
			out_file_name = target_name + "_vs_" + control_name + "_DE_elements_fc" + str(fc) + "_on_" + "_".join([reClass, reFamily, reName]) + ".txt"
			AnalyzeRNASeq.output_instances(de_members, summary, reClass, reFamily, reName, out_file_name)
		
			
			name = "_".join([reClass, reFamily, reName])
			#AnalyzeRNASeq.get_fc_histogram(summary, full_target_name, full_control_name, name, pc)
			AnalyzeRNASeq.get_log_fc_histogram(summary, full_target_name, full_control_name, name, target_library_size_vs_control_library_size, pc)
			# expression is in terms of full_control_name
			# fc is target/control
			AnalyzeRNASeq.get_expre_vs_fc(summary, full_target_name, full_control_name, name, target_library_size_vs_control_library_size, pc)
	
	
	output = open(target_name + "_derepressed_instances.pkl", 'wb')
	pickle.dump(de_members, output)
	output.close()
	
	output = open("RE_instance_counts.pkl", 'wb')
	pickle.dump(total_instances, output)
	output.close()
	
	#{reClass:{reFamily:reName:[id]}}}
	return (de_members, de_enrichment, total_instances)    
    
def find_K4up_and_K9down(de_members_K4, de_members_K9, total_instances, min_count=50):
	#de_members_K4:{reClass:{reFamily:reName:[id]}}}
	de_members_K4K9 = {}
	de_enrichment_K4K9 = {}
	for reClass in de_members_K4.keys():
		de_members_K4K9[reClass] = {}
		de_enrichment_K4K9[reClass] = {}
		for reFamily in de_members_K4[reClass].keys():
			de_members_K4K9[reClass][reFamily] = {}
			de_enrichment_K4K9[reClass][reFamily] = {}
			for reName in de_members_K4[reClass][reFamily].keys():
				K4up_ids = de_members_K4[reClass][reFamily][reName]
				K9down_ids = de_members_K9[reClass][reFamily][reName]
				intersection = list( set(K4up_ids) & set(K9down_ids) )
				total = total_instances[reClass][reFamily][reName]
				enrichment = len(intersection)*1.0/total
				de_members_K4K9[reClass][reFamily][reName] = intersection
				de_enrichment_K4K9[reClass][reFamily][reName] = enrichment
				
	
	mylist = AnalyzeRNASeq.rank_enrichment(de_members_K4K9, de_enrichment_K4K9, total_instances, "")
	out_file_name ="DE_enrichment_ranking_K4andK9_filtered_by_instancecount.txt"
	filter_ranked_list_by_RE_instance_count(mylist, total_instances, out_file_name, min_count)
	return (de_members_K4K9, de_enrichment_K4K9)

def find_K4up_or_K9down(de_members_K4, de_members_K9, total_instances, min_count=50):
	#de_members_K4:{reClass:{reFamily:reName:[id]}}}
	de_members_K4K9 = {}
	de_enrichment_K4K9 = {}
	for reClass in de_members_K4.keys():
		de_members_K4K9[reClass] = {}
		de_enrichment_K4K9[reClass] = {}
		for reFamily in de_members_K4[reClass].keys():
			de_members_K4K9[reClass][reFamily] = {}
			de_enrichment_K4K9[reClass][reFamily] = {}
			for reName in de_members_K4[reClass][reFamily].keys():
				K4up_ids = de_members_K4[reClass][reFamily][reName]
				K9down_ids = de_members_K9[reClass][reFamily][reName]
				myunion = list( set(K4up_ids).union(set(K9down_ids)) )
				total = total_instances[reClass][reFamily][reName]
				enrichment = len(myunion)*1.0/total
				de_members_K4K9[reClass][reFamily][reName] = myunion
				de_enrichment_K4K9[reClass][reFamily][reName] = enrichment
	
	mylist = AnalyzeRNASeq.rank_enrichment(de_members_K4K9, de_enrichment_K4K9, total_instances, "")
	out_file_name ="DE_enrichment_ranking_K4orK9_filtered_by_instancecount.txt"
	filter_ranked_list_by_RE_instance_count(mylist, total_instances, out_file_name, min_count)
	
	return (de_members_K4K9, de_enrichment_K4K9)	
	
def load_ids(file_name):
	id_list = []
	inf = open(file_name, "r")
	for line in inf:
		sline = line.strip()
		sline = sline.split()
		id_list.append(sline[0])
	return id_list

def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)


def analyze_subset(file_name, RLTR4mm):
	"""
	RLTR4MM: {id:
	{'annotation':RepElement class instance
	'H3K4me3_WT-W200-G200_rc':
	'H3K4me3_WT-W200-G200_rpkm':
	'H3K4me3_mir34bc_KO-W200-G200_rc':
	'H3K4me3_mir34bc_KO-W200-G200_rpkm':
	'H3K9me3_WT-W200-G400_rc':
	'H3K9me3_WT-W200-G400_rpkm':
	'H3K9me3_mir34bc_KO-W200-G400_rc':
	'H3K9me3_mir34bc_KO-W200-G400_rpkm':	
	}
}
	"""
	#print "\n load 5' RLTR4-mm"
	#file_name = "5_boundary_elements_on_2_LTR_derepressed_assembled_REs"
	
	mylist = load_ids(file_name)
	print "\n"
	print file_name
	print "Number of elements: ", len(mylist) 
	mydic = extract(RLTR4mm, mylist)
	
	feature_name = "H3K4me3_WT"
	full_feature_name_K4 =  feature_name + "-W200-G200_rpkm"
	H3K4me3_WT_enrichment = calculate_enrichment(mydic, full_feature_name_K4)
	print "WT H3K4me3 enrichment: ", H3K4me3_WT_enrichment
	
	feature_name = "H3K9me3_WT"
	full_feature_name_K9 =  feature_name + "-W200-G400_rpkm"
	H3K9me3_WT_enrichment = calculate_enrichment(mydic, full_feature_name_K9)
	print "WT H3K9me3 enrichment: ", H3K9me3_WT_enrichment
	
	present_list = [full_feature_name_K4, full_feature_name_K9]
	absent_list = []
	enrichment, positive_ids = calculate_pattern_enrichment_in_single_species(mydic, present_list, absent_list, threshold = 0.0001) 
	print "WT H3K4me3 and H3K9me3 co-enrichment: ", enrichment
	
	fc = 1.5
	min_rpkm = 0.01
	min_rc = 5
	pc = 5

	target_name = 'H3K4me3_mir34bc_KO'
	control_name =  "H3K4me3_WT"
	suffix = '-W200-G200'
	
	K4_target_library_size_vs_control_library_size = AnalyzeRNASeq.get_target_library_size_vs_control_library_size_from_single_RE_species (RLTR4mm, target_name + suffix, control_name + suffix)
	
	if K4_target_library_size_vs_control_library_size == -1:
		print "There is not a single locus with reads in both %s and %s" %(target_name, control_name)
		exit(1)
	
	name = "_on_" +  file_name
	AnalyzeRNASeq.get_log_fc_histogram(mydic, target_name + suffix, control_name + suffix, name, K4_target_library_size_vs_control_library_size, pc)
	AnalyzeRNASeq.get_expre_vs_fc(mydic,  target_name + suffix, control_name + suffix, name, K4_target_library_size_vs_control_library_size, pc)
	
	#enrichment_members: [id]
	#enrichment_statistics: ratio
	(enrichment_members_K4UP, enrichment_statistics, total) = AnalyzeRNASeq.get_upregulated_REs(mydic, target_name + suffix, control_name + suffix, fc, K4_target_library_size_vs_control_library_size, min_rpkm, min_rc, pc)
	print "Enrichment of copies of H3K4me3 up upon KO: ", enrichment_statistics, " of ", total, " copies"
	
	fc = 1.3
	min_rc = 3
	pc = 3

	target_name = "H3K9me3_WT" 
	control_name = 'H3K9me3_mir34bc_KO'
	suffix = '-W200-G400'
	K9_target_library_size_vs_control_library_size = AnalyzeRNASeq. get_target_library_size_vs_control_library_size_from_single_RE_species(RLTR4mm, target_name + suffix, control_name + suffix)
	
	if K9_target_library_size_vs_control_library_size == -1:
		print "There is not a single locus with reads in both %s and %s" %(target_name, control_name)
		exit(1)
	
	
	name = "_on_" +  file_name
	AnalyzeRNASeq.get_log_fc_histogram(mydic, target_name + suffix, control_name + suffix, name, K9_target_library_size_vs_control_library_size, pc)
	AnalyzeRNASeq.get_expre_vs_fc(mydic,  target_name + suffix, control_name + suffix, name, K9_target_library_size_vs_control_library_size, pc)
	
	(enrichment_members_K9DOWN, enrichment_statistics, total) = AnalyzeRNASeq.get_upregulated_REs(mydic, target_name + suffix, control_name + suffix, fc, K9_target_library_size_vs_control_library_size, min_rpkm, min_rc, pc)
	
	print "Enrichment of copies of H3K9me3 down upon KO: ", enrichment_statistics, " of ", total, " copies"
	
	
	intersection = list( set(enrichment_members_K4UP) & set(enrichment_members_K9DOWN) )
	enrichment = len(intersection)*1.0/total
	print "Intersection of K4 up and K9 down upon KO, enrichment: ", enrichment
	
	union = list( set(enrichment_members_K4UP).union (set(enrichment_members_K9DOWN)) )
	enrichment = len(union)*1.0/total
	print "Union:  K4 up or K9 down upon KO, enrichment: ", enrichment
	
	return enrichment_members_K4UP, enrichment_members_K9DOWN
    
def main(argv):
	parser = OptionParser()
	parser.add_option("-n", "--name_for_all_pickle_files", action="store", type="string", dest="summary_name", help="common name of all pickle files, one pickle for one reName", metavar="<str>")
	parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="RE tree in pickle format")
	parser.add_option("-l", "--RE_annotation_pkl_file_location", action="store", type="string", dest="RE_pkl_file_location", metavar="<str>", help="location of RE pkl files named in repClass_repFamily_repName.txt")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	#load the RE tree to get the RE file names
	print "\n\nLoading RE tree"
	re_tree = pickle.load(open(opt.RE_Tree, 'rb'))
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	
	print "\nLoading the RE instance counts"
	# {reClass:{reFamily:{reName:{feature_name:count}}}}
	RE_instance_counts = pickle.load(open("RE_instance_counts.pkl", 'rb') )
	
	# Test the redundancy in RE names
	print "\n\nTesting redundancy in RE names"
	name_set = set([])
	reused_name_set = set([])
	for reClass in re_tree.keys():
		for reFamily in re_tree[reClass].keys():
			overlap = name_set & set(re_tree[reClass][reFamily])
			reused_name_set = reused_name_set.union(overlap)
			name_set = name_set.union ( set(re_tree[reClass][reFamily]))
	name_list = list(name_set)
	reused_name_list = list(reused_name_set)
	print "There are (%d, %d) (singular,reused) names" %(len(name_list), len(reused_name_list))
	if (len(name_list) == numb_names):
		print "REs do not have redundant names"
	else:
		print "There is redundancy of names in REs:"
		print reused_name_list
	
	#load the annotation files and add annotation to features if that is not done
	print "\n\nExamine features and load annotation if necessary"
	index = 0
	for reClass in re_tree.keys():
		for reFamily in re_tree[reClass].keys():
			for reName in re_tree[reClass][reFamily]:
				summary_file_name = opt.summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				summary = pickle.load(inf)
				inf.close()
				index += 1
				feature_names = AnalyzeRNASeq.get_feature_names(summary) 
				if index <2: #only print once
					print "The features stored in the ", summary_file_name, " are ",  feature_names
				if 'annotation' not in feature_names:
					print "updating %s" %summary_file_name
					AnalyzeRNASeq.add_annotation_to_features(opt.RE_pkl_file_location, summary, opt.summary_name, reClass, reFamily, reName)
	
	current_dir = os.getcwd()
	
	print "\nGet target_library_size_vs_control_library_size"
	target_name = 'H3K4me3_mir34bc_KO'
	control_name =  "H3K4me3_WT"
	suffix = '-W200-G200'
	
	K4_target_library_size_vs_control_library_size = AnalyzeRNASeq.get_target_library_size_vs_control_library_size(re_tree, opt.summary_name, target_name + suffix, control_name + suffix)
	
	target_name = "H3K9me3_WT" 
	control_name = 'H3K9me3_mir34bc_KO'
	suffix = '-W200-G400'
	K9_target_library_size_vs_control_library_size = AnalyzeRNASeq.get_target_library_size_vs_control_library_size(re_tree, opt.summary_name, target_name + suffix, control_name + suffix)
	
	
	print "\nTest the hypothesis that the Histone mark levels in WT predicts the dereperssion"
	#feature_level: {reClass:{reFamily:{reName:{feature_name:level}}}}
	feature_level_mean, feature_level_median, feature_enrichment = get_feature_level(re_tree, opt.summary_name)
	
	#pickle_file_name ="summary_feature_level_mean.pkl"
	#output = open(pickle_file_name, 'wb')
	#pickle.dump(feature_level_mean, output)
	#output.close()
	
	#pickle_file_name ="summary_feature_level_median.pkl"
	#output = open(pickle_file_name, 'wb')
	#pickle.dump(feature_level_median, output)
	#output.close()
	
	pickle_file_name ="summary_feature_enrichment.pkl"
	output = open(pickle_file_name, 'wb')
	pickle.dump(feature_enrichment, output)
	output.close()
	
	#print "\nRank the RE species according to H3K4m3_WT level" 
	#feature_name = "H3K4me3_WT"
	#full_feature_name =  feature_name + "-W200-G200_rpkm"
	#out_file_name = feature_name + "_level_ranking.txt"
	#species_ranked_by_H3K4me3WT_level = rank_by_feature_level (feature_level_mean, full_feature_name, out_file_name)
	
	#print "\nRank the RE species according to H3K9m3_WT level" 
	#feature_name = "H3K9me3_WT"
	#full_feature_name =  feature_name + "-W200-G400_rpkm"
	#out_file_name = feature_name + "_level_ranking.txt"
	#species_ranked_by_H3K9me3WT_level = rank_by_feature_level (feature_level_mean, full_feature_name, out_file_name)
	
	print "\nOnly RE species with instance count >=50 are reported"
	min_count = 50
	
	feature_enrichment = pickle.load(open("summary_feature_enrichment.pkl", 'rb'))
	print "\nRank the RE species according to H3K4m3_WT enrichment" 
	feature_name = "H3K4me3_WT"
	full_feature_name =  feature_name + "-W200-G200_rpkm"
	out_file_name = feature_name + "_enrichment_ranking.txt"
	species_ranked_by_H3K4me3WT_enrichment = rank_by_feature_level (feature_enrichment, full_feature_name, "")
	filter_ranked_list_by_RE_instance_count(species_ranked_by_H3K4me3WT_enrichment, RE_instance_counts, out_file_name, min_count)
	
	print "\nRank the RE species according to H3K9m3_WT enrichment" 
	feature_name = "H3K9me3_WT"
	full_feature_name =  feature_name + "-W200-G400_rpkm"
	out_file_name = feature_name + "_enrichment_ranking.txt"
	species_ranked_by_H3K9me3WT_enrichment = rank_by_feature_level (feature_enrichment, full_feature_name, "")
	filter_ranked_list_by_RE_instance_count(species_ranked_by_H3K9me3WT_enrichment, RE_instance_counts, out_file_name, min_count)
	
	print "\nRank the RE species according to H3K4me3 and H3K9me3 co-enrichment"
	feature_name = "H3K4me3_WT"
	full_feature_name_K4 =  feature_name + "-W200-G200_rpkm"
	feature_name = "H3K9me3_WT"
	full_feature_name_K9 =  feature_name + "-W200-G400_rpkm"
	present_list = [full_feature_name_K4, full_feature_name_K9]
	absent_list = []
	#pattern_enrichment:{reClass:{reFamily:{reName:{feature_name:enrichment_ratio}}}}
	#pattern_positve_ids:{reClass:{reFamily:{reName:{feature_name:[ids]}}}}
	pattern_enrichment, pattern_positve_ids = find_pattern(re_tree, opt.summary_name, present_list, absent_list)
	species_ranked_by_K4K9WT_co_enrichment = rank_by_pattern_enrichment(pattern_enrichment, "") #[(reClass, reFamily, reName, enrichment)]
	out_file_name = "K4K9_WT_co" + "_enrichment_ranking.txt"
	min_count = 50
	filtered_species_ranked_by_K4K9WT_co_enrichment = filter_ranked_list_by_RE_instance_count(species_ranked_by_K4K9WT_co_enrichment, RE_instance_counts, out_file_name, min_count)
	
	#Output ids that exhibit the co-enrichment pattern of K4 and K9
	#reClass = "LINE"
	#reFamily = "L1"
	#reName = "L1Md_Gf"
	#out_filename = reClass + "_" + reFamily + "_" + reName + "_K4K9_ids"
	#output_pattern_ids(pattern_positve_ids, reClass, reFamily, reName, out_filename)
	
	
	
	fc = 1.5
	min_rpkm = 0.01
	min_rc = 5
	pc = 5
	print "\nCalculating H3K4me3 up in KO"
	target_name = 'H3K4me3_mir34bc_KO'
	control_name =  "H3K4me3_WT"
	suffix = '-W200-G200'
	(de_members_K4, de_enrichment_K4, total_instances) = find_differentially_marked_families(target_name, control_name, suffix, re_tree, opt.summary_name, fc, K4_target_library_size_vs_control_library_size, min_rpkm, min_rc, pc)
	
	fc = 1.5
	min_rc = 3
	pc = 3
	print ""
	print "Calculating H3K9me3 down in KO"
	target_name = "H3K9me3_WT" 
	control_name = 'H3K9me3_mir34bc_KO'
	suffix = '-W200-G400'
	(de_members_K9, de_enrichment_K9, total_instances) = find_differentially_marked_families(target_name, control_name, suffix, re_tree, opt.summary_name, fc, K9_target_library_size_vs_control_library_size, min_rpkm, min_rc, pc)
	
	print "\n intersect K4 and K9"
	find_K4up_and_K9down(de_members_K4, de_members_K9, total_instances)
	print "\n union K4 and K9"
	find_K4up_or_K9down(de_members_K4, de_members_K9, total_instances)
	
	print "\nLoad the RLTR4-mm summary"
	current_dir = os.getcwd()
	path = "/home/data/mm9/Lin/processed/Epigenome/summary"
	os.chdir(path)
	boundary_elements_file_name = "summary_on_LTR_ERV1_RLTR4_Mm.pkl"
	assert( Utility_extended.fileExists(boundary_elements_file_name) == 1)
	inf = open(boundary_elements_file_name, 'rb')
	RLTR4mm = pickle.load(inf)
	inf.close()
	print "There are %d %s as boundary elements" %(len(RLTR4mm.keys()), boundary_elements_file_name)
	
	
	os.chdir(current_dir)
	
	file_name = "5_boundary_elements_on_2_LTR_derepressed_assembled_REs"
	analyze_subset(file_name, RLTR4mm)
	
	file_name = "5_boundary_elements_on_2_LTR_nonderepressed_assembled_REs"
	analyze_subset(file_name, RLTR4mm)
	
	file_name = "3_boundary_elements_on_2_LTR_derepressed_assembled_REs"
	enrichment_members_K4UP, enrichment_members_K9DOWN = analyze_subset(file_name, RLTR4mm)
	print enrichment_members_K4UP
	
		
if __name__ == "__main__":
	main(sys.argv)
