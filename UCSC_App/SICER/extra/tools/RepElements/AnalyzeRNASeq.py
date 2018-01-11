#!/usr/bin/env python

"""
summary pickle data structure 
{id:
	{'annotation':RepElement class instance
	'index12_rc':
	'index12_rpkm':
	'index3_rc':
	'index3_rpkm':
	'index6_rc':
	'index6_rpkm':
	'index8_rc':
	'index8_rpkm':
	'index9_rc':
	'index9_rpkm':
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

def get_feature_names(summary):
	"""
	summary: {id:{feature_names: values}}
	return a list
	"""
	myids = summary.keys()
	myid = myids[0]
	feature_names = summary[myid].keys()
	feature_names.sort()
	return feature_names

def add_annotation_to_features(pkl_file_location, summary, summary_name, reClass, reFamily, reName):
	"""
	annotation file is in pkl_file_name = "_".join([reClass, reFamily, reName]) + ".pkl"
		data structure: a dictionary {id: rep_element class instance}
	summary:{id:{feature_names: values}}
	"""
	# load in all rep element instances from "_".join([reClass, reFamily, reName]) + ".pkl"
	currentdir = os.getcwd()
	os.chdir(pkl_file_location)
	pkl_file_name = "_".join([reClass, reFamily, reName]) + ".pkl"
	assert( Utility_extended.fileExists(pkl_file_name) == 1)
	inf = open(pkl_file_name, 'rb')
	MyRepElements = pickle.load(inf) # {id:Rep_Element class instance}
	inf.close()
	os.chdir(currentdir)
	
	# add the annotation into the summary
	for myid in summary.keys():
		summary[myid]['annotation'] = MyRepElements[myid]
	
	# write the updated summary back to 
	summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
	output = open(summary_file_name, 'wb')
	pickle.dump(summary, output)
	output.close()	

def get_subset(summaries, idset=None):
	"""
	summaries: [{id:{feature_names: values}}]
	return: [{id:{feature_names: values}}] with only ids in idset included.
	"""
	mysummaries = []
	if idset is not None:
		for summary in summaries:
			mykeys = list(set(summary.keys()) & set(idset))
			mysummary = {}
			for myid in mykeys:
				mysummary[myid] = copy.deepcopy(summary[myid])
			mysummaries.append(mysummary)
	else:
		mysummaries = copy.deepcopy(summaries)
	return mysummaries

def get_target_library_size_vs_control_library_size_from_single_RE_species(summary, target_name, control_name):
	target_library_size_vs_control_library_size = -1
	for myid in summary.keys():
		myre = summary[myid]
		if myre[target_name + "_rc"] > 0 and myre[control_name + "_rc"] > 0:
			rc_t = myre[target_name + "_rc"] * 1.0
			rc_c = myre[control_name + "_rc"] * 1.0
			rpkm_t = myre[target_name + "_rpkm"] * 1.0
			rpkm_c = myre[control_name + "_rpkm"] * 1.0
			target_library_size_vs_control_library_size = (rc_t/rc_c) / (rpkm_t/rpkm_c)
			print "The ratio between the library size of %s and the library size of %s is %f" %(target_name, control_name, target_library_size_vs_control_library_size)
			return target_library_size_vs_control_library_size
	return target_library_size_vs_control_library_size
	
def get_target_library_size_vs_control_library_size(re_tree, summary_name, target_name, control_name):
	"""
	We need to get the ratio of target_library_size over control_library_size
		rpkm = (rc /(library_size/1M) ) / (length/1k)
		rpkm_t/rpkm_c  = (rc_t/rc_c) / (library_size_t/library_size_c)
		(library_size_t/library_size_c) = (rc_t/rc_c) / (rpkm_t/rpkm_c)
		fc = (rcp_t/(library_size_t) / (rcp_c/library_size_c) ) = (rcp_t/rcp_c) / (library_size_t/library_size_c)
		
	If target_library_size_vs_control_library_size = -1 that means there is not a single locus that has non-zero in target and control at the same time  in summary
	"""
	target_library_size_vs_control_library_size = -1
	
	for reClass in re_tree.keys():
		for reFamily in re_tree[reClass].keys():
			for reName in re_tree[reClass][reFamily]:
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				summary = pickle.load(inf)
				inf.close()
				target_library_size_vs_control_library_size = get_target_library_size_vs_control_library_size_from_single_RE_species(summary, target_name, control_name)
				if target_library_size_vs_control_library_size > 0:
					return target_library_size_vs_control_library_size
	return target_library_size_vs_control_library_size
	
	
def get_expression_fc_for_RE_clusters(summaries, target_name, control_name, target_library_size_vs_control_library_size, pc = 1, idset=None):
	"""
	summaries: [{id:{feature_names: values}}]
	return {id:fc}
	"""
	
	mysummaries = get_subset(summaries, idset)
	
	fc = {}	#{id:fc}
	for mysummary in mysummaries:
		for myid in mysummary:
			myre = mysummary[myid]
			myfc = (myre[target_name + "_rc"] + pc)/float(myre[control_name + "_rc"] + pc)
			myfc = myfc / target_library_size_vs_control_library_size # Normalize by the library size ratio
			fc[myid] = myfc
	return fc

def get_rc_for_RE_cluster(summaries, target_name, control_name, idset=None):	
	"""
	summaries: [{id:{feature_names: values}}]
	return 
	(target_rc, control_rc)
	}
	"""
	mysummaries = get_subset(summaries, idset)
	
	target_rc = 0
	control_rc = 0
	for mysummary in mysummaries:
		for myid in mysummary.keys():
			myre = mysummary[myid]
			target_rc += myre[target_name + "_rc"]
			control_rc += myre[control_name + "_rc"]
	return (target_rc, control_rc)

def get_age_for_RE_cluster(summaries, mode, idset=None):
	"""
	summaries: [{id:{feature_names: values}}]
	"mode": "max", "min", "median", "mean"
	
	"""
	assert (mode=="max" or mode=="min" or mode == "mean" or mode == "median")
	
	mysummaries = get_subset(summaries, idset)
		
	age_list = []
	for mysummary in mysummaries:
		for myid in mysummary.keys():
			re_element = mysummary[myid]["annotation"]
			age_list.append(re_element.age)
		
	if mode == "max":
		my_age = max(age_list)
	elif mode == "min":
		my_age = min(age_list)
	elif mode == "mean":
		my_age = numpy.average(age_list)
	elif mode == "median":
		my_age = numpy.median(age_list)
	return my_age
	
def get_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1):
	"""
	"""
	fc_list = (get_expression_fc_for_RE_clusters([summary], target_name, control_name, target_library_size_vs_control_library_size, pc)).values()
	
	plt.hist(fc_list, bins=50, color='r', normed=True, alpha=0.75, label= target_name + "_vs_" + control_name + "_fc" + "_on_" + name )	
	plt.title(target_name + "_vs_" + control_name + "_on_" + name + "_fc_hist")
	plt.xlabel("fc")
	plt.ylabel("Frequency")
	plt.legend()
	#plt.legend(loc = 'upper left')
	plt.savefig(target_name + "_vs_" + control_name + "_on_" + name + "_fc_hist" + ".png", format="png")
	plt.close()
	return fc_list	

def get_log_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1):
	fc_list = (get_expression_fc_for_RE_clusters([summary], target_name, control_name, target_library_size_vs_control_library_size, pc)).values()
	log2_fc_list  = [log(item, 2) for item in fc_list]
	plt.hist(log2_fc_list, bins=50, color='r', normed=True, alpha=0.75, label= target_name + "_vs_" + control_name + "_fc" + "_on_" + name )	
	plt.title(target_name + "_vs_" + control_name + "_on_" + name + "_fc_hist")
	plt.xlabel("log(fc,2)")
	plt.ylabel("Frequency")
	plt.legend()
	#plt.legend(loc = 'upper left')
	plt.savefig(target_name + "_vs_" + control_name + "_on_" + name + "_logfc_hist" + ".png", format="png")
	plt.close()
	return log2_fc_list

	
def get_expre_vs_fc(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1):
	"""
	MA plot, x axis is expression rpkm in control, y axis is fc
	
	"""
	exp = []
	fc = []
	
	for myid in summary.keys():
		myre = summary[myid]
		myfc = (myre[target_name + "_rc"] + pc)/float(myre[control_name + "_rc"] + pc)
		myfc = myfc / target_library_size_vs_control_library_size
		
		exp.append(myre[control_name + "_rpkm"])
		fc.append(myfc)
	
	title = target_name + "_vs_" + control_name + "_on_" + name + "_exp_vs_fc"
	xscale = "log"
	yscale = "log"
	plt.plot(exp, fc, "bo", markersize = 1.5, alpha = 1);
	plt.xlabel("exp")
	plt.ylabel("fc")
	ax = plt.gca();
	ax.set_xscale(xscale)
	ax.set_yscale(yscale)
	#ax.set_aspect(1.) 
	ax.grid (color='gray', linestyle='dashed')
	#plt.savefig(title + ".eps", format="eps")
	plt.savefig(title + ".png", format="png")
	plt.close()
	
def get_upregulated_REs(summary, target_name, control_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc=15, pc = 1):
	"""
	Find re instances with fold change beyone fc
	
	target_name: eg, index8 or index6
	summary: {id:{feature:value}}
	
	previously: 
	fc refers to target_name + "_rc" vs control_name + "_rc"
	this is wrong because the library size might be different.
	
	We need to get the ratio of target_library_size over control_library_size
		rpkm = (rc /(library_size/1M) ) / (length/1k)
		rpkm_t/rpkm_c  = (rc_t/rc_c) / (library_size_t/library_size_c)
		(library_size_t/library_size_c) = (rc_t/rc_c) / (rpkm_t/rpkm_c)
		fc = (rcp_t/(library_size_t) / (rcp_c/library_size_c) ) = (rcp_t/rcp_c) / (library_size_t/library_size_c)
	
	enrichment_members: [id]
	enrichment_statistics: ratio
	total: total_number of elements 
	"""
	enrichment_members = []
	enrichment_statistics = 0.0
	total = len(summary.keys())
	num_de = 0
	
	for myid in summary.keys():
		myre = summary[myid]
		if myre[target_name + "_rc"] >= min_rc and myre[target_name + "_rpkm"] >= min_rpkm:
			myfc = (myre[target_name + "_rc"] + pc)/float(myre[control_name + "_rc"] + pc)
			myfc = myfc / target_library_size_vs_control_library_size
			if myfc >= fc:
				num_de += 1
				enrichment_members.append(myid)
	enrichment_statistics = float(num_de)/float(total)
	return (enrichment_members, enrichment_statistics, total)	

def get_all_upregulated_REs(re_tree, summary_name, target_name, control_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc=15, pc = 1):
	"""
	
	"""
	de_members = {} #{reClass:{reFamily:reName:[id]}}}
	de_enrichment = {}  #{reClass:{reFamily:reName:enrichment_ratio}}}
	total_instances = {} #{reClass:{reFamily:{reName:count}}}
	for reClass in re_tree.keys():
		de_members[reClass] = {}
		de_enrichment[reClass] = {}
		total_instances[reClass] = {}
		for reFamily in re_tree[reClass].keys():
			de_members[reClass][reFamily] = {}
			de_enrichment[reClass][reFamily] = {}
			total_instances[reClass][reFamily] = {}
			for reName in re_tree[reClass][reFamily]:
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				summary = pickle.load(inf)
				inf.close()
				#de_re_instances:[id]
				#de_re_enrichment: ratio
				#total:total_number of elements 
				(de_re_instances, de_re_enrichment, total) = get_upregulated_REs(summary, target_name, control_name, fc, target_library_size_vs_control_library_size,  min_rpkm, min_rc, pc)
				de_members[reClass][reFamily][reName] = de_re_instances
				de_enrichment[reClass][reFamily][reName] = de_re_enrichment
				total_instances[reClass][reFamily][reName] = total
	return (de_members, de_enrichment, total_instances)

def family_level_number_of_instances(total_instances):
	"""
	return: family_level_total_counts[reClass][reFamily]: number of instances in this family
	"""
	family_level_total_counts = {}
	for reClass in total_instances.keys():
		family_level_total_counts[reClass] = {}
		for reFamily in total_instances[reClass].keys():
			family_level_total_counts[reClass][reFamily] = 0
			total = 0
			for reName in total_instances[reClass][reFamily].keys():
				total += total_instances[reClass][reFamily][reName]
			family_level_total_counts[reClass][reFamily] = total
	return 	family_level_total_counts	

	
def get_name_level_enrichment_counts(de_members, reClass, reFamily):
	"""
	name_level_enrichment_counts[reFamily]: # of instances shown derepression
	"""
	name_level_enrichment_counts = {}
	for reName in de_members[reClass][reFamily].keys():
		name_level_enrichment_counts[reName]= len(de_members[reClass][reFamily][reName])
	return 	name_level_enrichment_counts

def get_entries_above_cutoff(dic, cutoff=0):
	above_cutoff_dic = {}
	for myid in dic.keys():
		if dic[myid] > cutoff:
			above_cutoff_dic[myid] = dic[myid]
	return above_cutoff_dic
	
def get_name_level_enrichment_ratios(de_members, reClass, reFamily, total_instances):
	"""
	If cutoff is 0, pick non-zero items
	"""
	name_level_enrichment_counts = get_name_level_enrichment_counts(de_members, reClass, reFamily)
	name_level_enrichment_ratios = {}
	for reName in name_level_enrichment_counts.keys():
		if total_instances[reClass][reFamily][reName] == 0:
			name_level_enrichment_ratios[reName] = 0
		else:
			name_level_enrichment_ratios[reName] = name_level_enrichment_counts[reName]/float(total_instances[reClass][reFamily][reName])
	return name_level_enrichment_ratios
		
def get_family_level_enrichment_counts(de_members, reClass):
	"""
	family_level_enrichment_counts[reFamily]: # of instances shown derepression
	"""
	family_level_enrichment_counts = {}
	for reFamily in de_members[reClass].keys():
		family_level_enrichment_counts[reFamily] = 0
		total = 0
		for reName in de_members[reClass][reFamily].keys():
				total += len(de_members[reClass][reFamily][reName])
		family_level_enrichment_counts[reFamily] = total
	return 	family_level_enrichment_counts

def get_family_level_enrichment_ratios(de_members, reClass, total_instances, cutoff = 0):
	"""
	de_members:  #{reClass:{reFamily:reName:[id]}}}
	total_instances: #{reClass:{reFamily:reName:number of instances}}}
	family_level_total_counts: from family_level_number_of_instances
	family_level_enrichment_ratios[reFamily]: ratio of derepressed instances
	"""
	family_level_enrichment_counts= get_family_level_enrichment_counts (de_members, reClass)
	family_level_total_counts = family_level_number_of_instances(total_instances)
	family_level_enrichment_ratios = {}
	for reFamily in family_level_enrichment_counts.keys():
		total_number_instances = family_level_total_counts[reClass][reFamily]
		if  total_number_instances < 0.01:
			family_level_enrichment_ratios[reFamily] = 0
		else:
			family_level_enrichment_ratios[reFamily] = float(family_level_enrichment_counts[reFamily])/total_number_instances
		
	above_threshold_family_level_enrichment_ratios = {}
	for reFamily in family_level_enrichment_ratios.keys():
		if family_level_enrichment_counts[reFamily] >= cutoff:
			above_threshold_family_level_enrichment_ratios[reFamily] = family_level_enrichment_ratios[reFamily]
	return 	family_level_enrichment_ratios, above_threshold_family_level_enrichment_ratios
	
def get_class_level_number_of_instances(total_instances):
	"""
	class_level_total_counts[reClass]: number of instances in this family
	"""
	class_level_total_counts = {}
	for reClass in total_instances.keys():
		class_level_total_counts[reClass] = 0
		total = 0
		for reFamily in total_instances[reClass].keys():
			for reName in total_instances[reClass][reFamily].keys():
				total += total_instances[reClass][reFamily][reName]
		class_level_total_counts[reClass] = total
	return 	class_level_total_counts	
	
def get_class_level_enrichment_counts(de_members):
	"""
	de_members:{reClass:{reFamily:reName:[id]}}}
	class_level_enrichment_counts[reClass]: the number of derepressed instances in the class
	"""
	class_level_enrichment_counts = {}
	for reClass in de_members.keys():
		class_level_enrichment_counts[reClass] = 0
		total = 0
		for reFamily in de_members[reClass].keys():
			for reName in de_members[reClass][reFamily].keys():
				total += len(de_members[reClass][reFamily][reName])
		class_level_enrichment_counts[reClass] = total
	
	return class_level_enrichment_counts
				
def get_class_level_enrichment_ratios(de_members, total_instances, cutoff=10):
	"""
	class_level_enrichment_counts[reClass]: the number of derepressed instances in the class
	non_zero_class_level_enrichment_counts: only record the non-zero entries, 
	cutoff: the minimum number of derepressed instances for inclusive in  
	"""
	class_level_enrichment_counts = get_class_level_enrichment_counts(de_members)
	class_level_number_of_instances = get_class_level_number_of_instances(total_instances)
	
	class_level_enrichment_ratios = {}
	for reClass in class_level_enrichment_counts.keys():
		if class_level_number_of_instances[reClass] <0.01:
			class_level_enrichment_ratios[reClass] = 0
		else:
			class_level_enrichment_ratios[reClass] = float(class_level_enrichment_counts[reClass])/class_level_number_of_instances[reClass]
	
	# only non-zero entries
	above_cutoff_class_level_enrichment_ratios = {}
	for reClass in class_level_enrichment_counts.keys():
		if class_level_enrichment_counts[reClass] > 0.0 and class_level_enrichment_counts[reClass] >= cutoff:
			above_cutoff_class_level_enrichment_ratios[reClass]= class_level_enrichment_ratios[reClass]
	return class_level_enrichment_ratios, above_cutoff_class_level_enrichment_ratios


	
	
def rank_enrichment(de_members, de_enrichment, total_instances, out_file_name =""):
	"""
	out_file_name = target_name + "_vs_" + control_name + "_DE_enrichment_fc" + str(fc) + "_on_REs.txt"
	mylist:[(reClass, reFamily, reName, de_enrichment[reClass][reFamily][reName], total_instances[reClass][reFamily][reName])]
	"""
	mylist = []
	for reClass in de_members.keys():
		for reFamily in de_members[reClass].keys():
			for reName in de_members[reClass][reFamily].keys():
				mylist.append((reClass, reFamily, reName, de_enrichment[reClass][reFamily][reName], total_instances[reClass][reFamily][reName]))
	mylist.sort(key = itemgetter(3), reverse=True)
	if out_file_name != "":
		outf = open(out_file_name, "w")
		for item in mylist:
			outf.write("\t".join(map(str,item)) + "\n")
		outf.close()
	return mylist
	
def output_instances(de_members, summary, reClass, reFamily, reName, out_file_name):
	"""
	de_members: #{reClass:{reFamily:reName:[id]}}}
	summary: the summary for a particular reClass, reFamily, reName
	"""
	feature_names = get_feature_names(summary)
	feature_names.remove('annotation')
	
	myids = de_members[reClass][reFamily][reName]
	outf = open(out_file_name, "w")
	outline = "\t".join(["id", "chrom", "strand", "start", "end"]) + "\t"
	outline  += "\t".join(feature_names) + "\n"
	outf.write(outline)
	for myid in myids:
		re_element = summary[myid]['annotation']
		chrom = re_element.chrom
		start = re_element.genoStart
		end = re_element.genoEnd
		strand =  re_element.strand
		outline = "\t".join([myid, chrom, strand, str(start), str(end)]) + "\t"
		outline  += "\t".join([str(summary[myid][feature_name]) for feature_name in feature_names]) + "\n"
		outf.write(outline)
	outf.close()

def output_all_instances(de_members, summary_name, out_file_name):
	"""
	de_members: #{reClass:{reFamily:reName:[id]}}}
	
	output info in all ids defined in de_members
	"""
	# load in all rep element instances 
	
	outf = open(out_file_name, "w")
	
	for reClass in de_members.keys():
		for reFamily in de_members[reClass].keys():
			for reName in de_members[reClass][reFamily]:
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				summary = pickle.load(inf)
				inf.close()
				feature_names = get_feature_names(summary) 
				feature_names.remove('annotation')
	
				myids = de_members[reClass][reFamily][reName]
				for myid in myids:
					re_element = summary[myid]['annotation']
					chrom = re_element.chrom
					start = re_element.genoStart
					end = re_element.genoEnd
					strand =  re_element.strand
					outline = "\t".join([myid, chrom, strand, str(start), str(end)]) + "\t"
					outline  += "\t".join([str(summary[myid][feature_name]) for feature_name in feature_names]) + "\n"
					outf.write(outline)
	outf.close()	

def find_mappable_REs(summary):
	"""
	Find number of reads on RE across all libraries
	
	return {id:rc}
	"""
	mykeys = (summary.keys())
	
	# find REs that are mappable. if total_rc[mykey] > 0, it is mappable
	features = ["index12_rc", "index6_rc", "index8_rc", "index9_rc", "index3_rc"]
	total_rc = {} #{id:rc}
	for mykey in mykeys:
		total_rc[mykey] = 0
	for mykey in mykeys:
		for feature_name in features:
			total_rc[mykey] += summary[mykey][feature_name]
	return total_rc
	
def age_vs_derepression_for_RE_family(summary, target_name, control_name, target_library_size_vs_control_library_size, rc_threshold = 3, pc = 1):
	"""
	summary: {id:{feature_name:value}}
	rc_threshold: the minimum amount of reads on a RE for that RE to be deemed mappable.
	Return: {id:(age, quality, myfc)}
	"""
	print get_feature_names(summary)
	mykeys = sorted(summary.keys())
	
	# find REs that are mappable. if total_rc[mykey] >= rc_threshold, it is mappable
	total_rc = find_mappable_REs(summary) #{id:rc}
	
	myresult = {}
	for myid in mykeys:
		myre = summary[myid]
		re_element = myre['annotation']
		age = re_element.age
		quality = re_element.quality
		myfc = (float(myre[target_name + "_rc"]) + pc)/(float(myre[control_name + "_rc"]) + pc)
		myfc = myfc/target_library_size_vs_control_library_size
		myresult[myid] = (age, quality, myfc)
		
	myresult_on_mappable_REs = {}
	for myid in mykeys:
		if total_rc[myid] >= rc_threshold:
			myresult_on_mappable_REs[myid] = myresult[myid]
	
	return (myresult, myresult_on_mappable_REs)

	
	
def age_vs_derepression(summary_name, target_name, control_name, selected_members, target_library_size_vs_control_library_size, rc_threshold = 3, pc = 1):
	"""
	selected_members: {reClass:{reFamily:reName:[id]}}}
	summary_name: the macro name for all pickle
	
	Return: {id:(age, quality, myfc)}
	"""
	
	myresult ={}
	
	for reClass in selected_members.keys():
		for reFamily in selected_members[reClass].keys():
			for reName in selected_members[reClass][reFamily]:
				summary_file_name = summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				assert( Utility_extended.fileExists(summary_file_name) == 1)
				inf = open(summary_file_name, 'rb')
				summary = pickle.load(inf)
				inf.close()
				
				myids = selected_members[reClass][reFamily][reName]
				for myid in myids:
					myre = summary[myid]
					re_element = myre['annotation']
					age = re_element.age
					quality = re_element.quality
					myfc = (float(myre[target_name + "_rc"]) + pc)/(float(myre[control_name + "_rc"]) + pc)
					myfc = myfc/target_library_size_vs_control_library_size
					myresult[myid] = (age, quality, myfc)
				
				# find REs that are mappable. if total_rc[mykey] > 0, it is mappable
				total_rc = find_mappable_REs(summary)#{id:rc}
				
				myresult_on_mappable_REs = {}
				for myid in mykeys:
					if total_rc[myid] >= rc_threshold:
						myresult_on_mappable_REs[myid] = myresult[myid]
				
	return (myresult, myresult_on_mappable_REs)



def explore_sandwich_structure_for_a_region (region, extension, boundary_elements_by_chrom, shift=0.01):
	"""
	region: (start, end)
	boundary_elements_by_chrom: {id:{feature_name:value}},  on the same chrom as re_element
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE. As a result: extension = min ( extension, length_of_RE -10)
	
	return the number of boundary elements present, + shift is for plotting visibility. 
	"""
	start = region[0]
	end = region[1]
	length = end - start + 1
	extension = max(0, min(length - 10, extension))
	region_list = [(boundary_elements_by_chrom[myid]["annotation"].genoStart -extension, boundary_elements_by_chrom[myid]["annotation"].genoEnd + extension) for myid in boundary_elements_by_chrom.keys() ]
	left_overlap = Utility_extended.is_inside(start, region_list)
	right_overlap = Utility_extended.is_inside(end, region_list)
	return left_overlap + right_overlap + shift #shift is for plotting visibility	

def explore_sandwich_structure_for_a_region_w_trace (region, extension, boundary_elements_by_chrom, boundary_element_clustering_extension = 20):
	"""
	region: (start, end)
	boundary_elements_by_chrom: {id:{feature_name:value}},  on the same chrom as re_element
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE. As a result: extension = min ( extension, length_of_RE -10)
	2*boundary_element_clustering_extension:  the max distance between two bundary elements so that they are treated as one cluster
	
	return left_boundary_element_cluster [ids]; right_boundary_element_cluster[ids]. 
	"""
	start = region[0]
	end = region[1]
	length = end - start + 1
	extension = max(0, min(length - 10, extension))
	region_list = [(boundary_elements_by_chrom[myid]["annotation"].genoStart -extension, boundary_elements_by_chrom[myid]["annotation"].genoEnd + extension, myid) for myid in boundary_elements_by_chrom.keys() ] # [(start, end, id)]
	
	temp = Utility_extended.find_covering_cluster_of_regions(start, region_list, boundary_element_clustering_extension) #[(start, end, id)]
	left_boundary_element_cluster = [item[2] for item in temp] #[id]
	temp = []
	temp = Utility_extended.find_covering_cluster_of_regions(end, region_list, boundary_element_clustering_extension)#[(start, end, id)]
	right_boundary_element_cluster =  [item[2] for item in temp]  #[ids]
	
	return left_boundary_element_cluster, right_boundary_element_cluster 
	
	
def explore_sandwich_structure_for_one_RE(re_element, extension, boundary_elements_by_chrom, shift = 0.01):
	"""
	re_element: a single re_element class object
	boundary_elements_by_chrom: {id:{feature_name:value}},  on the same chrom as re_element
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE. As a result: extension = min ( extension, length_of_RE -10)
	
	return the number of boundary elements present, + 0.01 is for plotting. 
	"""
	region = (re_element.genoStart, re_element.genoEnd)
	return explore_sandwich_structure_for_a_region (region, extension, boundary_elements_by_chrom, shift)

def explore_sandwich_structure_for_RE_family(summary, target_name, control_name, extension, boundary_elements, target_library_size_vs_control_library_size, rc_threshold = 3, pc = 1):
	"""
	test the sandwich structure of the REs for derepressed REs. The idea is that only complete REs are able to be derepressed
	Might not be the best place to put it here
	selected_members: #{reClass:{reFamily:reName:[id]}}}
	boundary_elements: {id:{feature_name:value}}
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE.As a result: extension = min ( extension, length_of_RE -10)
	
	Return: {id:(number of boundary elements, fc)}
			{id:(number of boundary elements, fc)}
	"""
	
	boundary_elements_by_chrom = separate_by_chrom(boundary_elements)
	
	mykeys = summary.keys()
	
	myresult = {}
	for myid in mykeys:
		myre = summary[myid]
		re_element = myre['annotation']
		chrom = re_element.chrom
		myfc = (float(myre[target_name + "_rc"]) + pc)/(float(myre[control_name + "_rc"]) + pc)
		myfc = myfc/target_library_size_vs_control_library_size
		if chrom not in boundary_elements_by_chrom.keys():
			myresult[myid] = (0, myfc)
		else:
			eligibility = explore_sandwich_structure_for_one_RE(re_element, extension, boundary_elements_by_chrom[chrom]) # colocalization structure test
			myresult[myid] = (eligibility, myfc)
			
	myresult_on_mappable_REs = {}
	total_rc = find_mappable_REs(summary)
	for myid in mykeys:
		if total_rc[myid] >= rc_threshold:
			myresult_on_mappable_REs[myid] = myresult[myid]
	
	return (myresult, myresult_on_mappable_REs)	

def explore_sandwich_structure_for_clusteredRE_family(summary, target_name, control_name,	cluster_extension,  extension, boundary_elements, target_library_size_vs_control_library_size, pc = 1):
	"""
	test the sandwich structure of the REs for derepressed REs. The idea is that only complete REs are able to be derepressed
	Might not be the best place to put it here
	cluster_extension: the extension of each RE on each side for clustering of REs
	boundary_elements: {id:{feature_name:value}}
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE.As a result: extension = min ( extension, length_of_RE -10)
	pc: pseudo count for fold change
	
	Output: 
	
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
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	
	
	"""
	return explore_sandwich_structure_for_clusteredRE_families([summary], target_name, control_name, cluster_extension, extension, boundary_elements, target_library_size_vs_control_library_size, pc)
	
	
def explore_sandwich_structure_for_clusteredRE_families(summaries, target_name, control_name, cluster_extension, extension, boundary_elements, target_library_size_vs_control_library_size, pc = 1):
	"""
	test the sandwich structure of the REs for derepressed REs. The idea is that only REs with are able to be derepressed
	summaries: [{id:{feature_name:value}}]
	cluster_extension: the extension of each RE on each side for clustering of REs
	boundary_elements: {id:{feature_name:value}}
	extension: the max distance for between the boundary element and the edge of the element in consideration. There is the issue that the element in consideration can be really short, in which case the extension must be shorter, otherwise a boundary element will be counted as neighbouring both sides of the RE.As a result: extension = min ( extension, length_of_RE -10)
	
	pc: pseudo count for fold change calculation
	
	Output: 
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
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	
	"""
	
	
	boundary_elements_by_chrom = separate_by_chrom(boundary_elements) #{chrom:{id:{feature_name:value}}}
	
	number_of_families = len(summaries)
	#print "number_of_families: ", number_of_families
	

	#{chrom: {(region_start, region_end):[(start, end, id)]}}
	assembled_REs = copy.deepcopy(Union_REs_in_families(summaries, cluster_extension))
	
	#Find boundary elements for each cluster
	total_num_clusters = 0
	num_p_strand_clusters = 0
	num_n_strand_clusters = 0
	num_I_strand_clusters = 0
	boundary_element_clustering_extension = 20
	shift = 0.01
	for chrom in assembled_REs.keys():
		total_num_clusters += len(assembled_REs[chrom])
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			num_boundary_elements = 0 
			if chrom not in boundary_elements_by_chrom.keys() or len(boundary_elements_by_chrom[chrom]) == 0:
				left_boundary_element_cluster = []
				right_boundary_element_cluster = []	
			else:
				left_boundary_element_cluster, right_boundary_element_cluster = explore_sandwich_structure_for_a_region_w_trace (mycluster, extension, boundary_elements_by_chrom[chrom], boundary_element_clustering_extension)
				if len(left_boundary_element_cluster) > 0 :
					num_boundary_elements += 1
				if len(right_boundary_element_cluster) > 0 :
					num_boundary_elements += 1		
			assembled_REs[chrom][mycluster]["num_boundary_elements"] = num_boundary_elements + shift
			
			if assembled_REs[chrom][mycluster]["strand"] == "+":
				assembled_REs[chrom][mycluster]["5_boundary_elements"] = left_boundary_element_cluster
				assembled_REs[chrom][mycluster]["3_boundary_elements"] = right_boundary_element_cluster
				assembled_REs[chrom][mycluster]["I_boundary_elements"] = []
				num_p_strand_clusters += 1
			elif assembled_REs[chrom][mycluster]["strand"] == "-":
				assembled_REs[chrom][mycluster]["5_boundary_elements"] = right_boundary_element_cluster
				assembled_REs[chrom][mycluster]["3_boundary_elements"] = left_boundary_element_cluster
				assembled_REs[chrom][mycluster]["I_boundary_elements"] = []
				num_n_strand_clusters += 1
			else: # assembled_REs[chrom][mycluster]["strand"] == "I":
				assembled_REs[chrom][mycluster]["5_boundary_elements"] = []
				assembled_REs[chrom][mycluster]["3_boundary_elements"] = []
				assembled_REs[chrom][mycluster]["I_boundary_elements"] = left_boundary_element_cluster + right_boundary_element_cluster
				num_I_strand_clusters += 1
			
			#add age to assembled_REs, use the extreme value to represent the cluster
			mode = "min"
			ids_in_mycluster = assembled_REs[chrom][mycluster]["elements"]
			age = get_age_for_RE_cluster(summaries, mode, ids_in_mycluster)
			assembled_REs[chrom][mycluster]["age"] = age
			
			five_boundary_element_ids = assembled_REs[chrom][mycluster]["5_boundary_elements"]
			if len(five_boundary_element_ids) > 0:
				five_boundary_elements_age = get_age_for_RE_cluster([boundary_elements], mode, five_boundary_element_ids)
			else:
				five_boundary_elements_age = 0
			assembled_REs[chrom][mycluster]["5_boundary_elements_age"] = five_boundary_elements_age
			
			three_boundary_element_ids = assembled_REs[chrom][mycluster]["3_boundary_elements"]
			if len(three_boundary_element_ids) > 0:
				three_boundary_elements_age = get_age_for_RE_cluster([boundary_elements], mode, three_boundary_element_ids)
			else:
				three_boundary_elements_age = 0
			assembled_REs[chrom][mycluster]["3_boundary_elements_age"] = three_boundary_elements_age
			
			#add expression_fc to assembled_REs, use the extreme value to represent the cluster
			expression_fc_for_ids_in_mycluster = get_expression_fc_for_RE_clusters(summaries, target_name, control_name, target_library_size_vs_control_library_size, pc, ids_in_mycluster) #{id:fc}
			# Use the most extreme change among members of the cluster as representative
			my_max = max(expression_fc_for_ids_in_mycluster.values())
			my_min = min(expression_fc_for_ids_in_mycluster.values())
			if abs(log(my_max)) >= abs(log(my_min)):
				expression_fc_extreme = my_max
			else:
				expression_fc_extreme = my_min
			assembled_REs[chrom][mycluster]["expression_fc_" + target_name + "_vs_" + control_name] = expression_fc_extreme
			
			#add rc 
			(target_rc, control_rc) = get_rc_for_RE_cluster(summaries, target_name, control_name, ids_in_mycluster)
			
			#print ids_in_mycluster, (target_rc, control_rc)
			assembled_REs[chrom][mycluster][target_name + "_rc"] = target_rc
			assembled_REs[chrom][mycluster][control_name + "_rc"] = control_rc
			
	print "There are %d clusters " %(total_num_clusters)
	print "%d clusters have strand + " %num_p_strand_clusters
	print "%d clusters have strand - " %num_n_strand_clusters
	print "%d clusters have strand I " %num_I_strand_clusters
	
	# save assembled_REs to pkl
	print "Saving assembled REs to pkl"
	print 
	assembled_REs_file_name = "Assembled_REs.pkl"
	output = open(assembled_REs_file_name, 'wb')
	pickle.dump(assembled_REs, output)
	output.close()
	
	return assembled_REs
	
def output_assembled_REs(assembled_REs, target_name, control_name, output_filename):
	"""
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
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	
	"""
	assembled_RE_list = []
	for chrom in assembled_REs.keys():
		for region in assembled_REs[chrom].keys():
			start = region[0]
			end = region[1]
			strand = assembled_REs[chrom][region]["strand"]
			num_boundary_elements = assembled_REs[chrom][region]["num_boundary_elements"]
			five_boundary_elements = ",".join([str(item) for item in assembled_REs[chrom][region]["5_boundary_elements"]])
			five_boundary_elements_age = assembled_REs[chrom][region]["5_boundary_elements_age"]
			three_boundary_elements = ",".join([str(item) for item in assembled_REs[chrom][region]["3_boundary_elements"]])
			three_boundary_elements_age = assembled_REs[chrom][region]["3_boundary_elements_age"]
			i_boundary_elements = ",".join([str(item) for item in assembled_REs[chrom][region]["I_boundary_elements"]])
			elements = ",".join(assembled_REs[chrom][region]["elements"])
			age = assembled_REs[chrom][region]["age"]
			expression_fc = assembled_REs[chrom][region]["expression_fc_"+target_name+"_vs_" + control_name]
			target_rc = assembled_REs[chrom][region][target_name + "_rc"]
			control_rc = assembled_REs[chrom][region][control_name + "_rc"]
			assembled_RE_list.append((start, end, chrom, strand, expression_fc, num_boundary_elements, five_boundary_elements, elements, three_boundary_elements, five_boundary_elements_age, age, three_boundary_elements_age, target_rc, control_rc, i_boundary_elements))
	assembled_RE_list.sort( key = itemgetter(5, 4), reverse=True)
	outf = open(output_filename, "w")
	outline = "#start\tend\tchrom\tstrand\texpression_fc\tnum_boundary_elements\t5_LTR\tint\t3_LTR\t5_LTR_age\tint_age\t3_LTR_age\t" + target_name  + "_rc\t" + control_name + "_rc\tIndeterminateBoundaryElements\n"
	outf.write(outline)
	for item in assembled_RE_list:
		str_item = [str(c) for c in item]
		outline = "\t".join(str_item) +"\n"
		outf.write(outline)
	outf.close()

def get_boundary_elements_from_assembled_REs(assembled_REs):
	"""
	return classified boundary elements [id]
	"""
	# Collect boundary element ids on boundaries vs not_on_boundaries
	
	total_num_clusters = 0
	all_5_boundary_elements=[]
	all_3_boundary_elements=[]
	all_I_boundary_elements=[]
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			if assembled_REs[chrom][mycluster]["5_boundary_elements"]:
				total_num_clusters += 1
			all_5_boundary_elements += assembled_REs[chrom][mycluster]["5_boundary_elements"]
			if assembled_REs[chrom][mycluster]["3_boundary_elements"]:
				total_num_clusters += 1
			all_3_boundary_elements += assembled_REs[chrom][mycluster]["3_boundary_elements"]
			all_I_boundary_elements += assembled_REs[chrom][mycluster]["I_boundary_elements"]
	print "There are %d boundary clusters, excluding indeterminate ones" %total_num_clusters
	print "There are (%d, %d, %d) (5, 3, indeterminate) boundary elements. Totaling %d" %(len(all_5_boundary_elements), len(all_3_boundary_elements), len(all_I_boundary_elements), len(all_5_boundary_elements) + len(all_3_boundary_elements) + len(all_I_boundary_elements))

	return all_5_boundary_elements, all_3_boundary_elements, all_I_boundary_elements

def get_fc_from_assembled_REs(assembled_REs, target_name, control_name):
	"""
	Return [fc]
	"""
	fc_list = []
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			fc_list.append(assembled_REs[chrom][mycluster]["expression_fc_"+target_name+"_vs_" + control_name])
	return fc_list

def get_age_from_assembled_REs(assembled_REs):
	"""
	Return [age]
	"""
	age_list = []
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			age_list.append(assembled_REs[chrom][mycluster]["age"])
	return age_list
	
def get_assembled_REs_subset_by_num_LTR(assembled_REs, num_LTR):
	my_assembled_REs = {}
	for chrom in assembled_REs.keys():
		my_assembled_REs[chrom] = {}
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			my_num_LTR = assembled_REs[chrom][mycluster]["num_boundary_elements"]
			if  my_num_LTR >= num_LTR and my_num_LTR < (num_LTR+1):
				my_assembled_REs[chrom][mycluster] = copy.deepcopy(assembled_REs[chrom][mycluster])
	return my_assembled_REs
	
def get_assembled_REs_subset_by_rc(assembled_REs, target_name, control_name, rc_threshold):
	my_assembled_REs = {}
	for chrom in assembled_REs.keys():
		my_assembled_REs[chrom] = {}
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			target_rc = assembled_REs[chrom][mycluster][target_name + "_rc"]
			control_rc = assembled_REs[chrom][mycluster][control_name + "_rc"]
			if  target_rc >= rc_threshold or control_rc >= rc_threshold:
				my_assembled_REs[chrom][mycluster] = copy.deepcopy(assembled_REs[chrom][mycluster])
	return my_assembled_REs


def classify_assembled_REs_by_fc(assembled_REs, target_name, control_name, fc_threshold):
	"""
	return differentially expressed, and non-differentially expressed
	"""
	my_assembled_REs = {}
	out_assembled_REs = {}
	for chrom in assembled_REs.keys():
		my_assembled_REs[chrom] = {}
		out_assembled_REs[chrom] = {}
		for mycluster in assembled_REs[chrom].keys(): #mycluster: (start,end)
			fc = assembled_REs[chrom][mycluster]["expression_fc_"+target_name+"_vs_" + control_name]
			if  fc >= fc_threshold:
				my_assembled_REs[chrom][mycluster] = copy.deepcopy(assembled_REs[chrom][mycluster])
			else:
				out_assembled_REs[chrom][mycluster] = copy.deepcopy(assembled_REs[chrom][mycluster])
	return my_assembled_REs, out_assembled_REs
	
	
def get_number_of_assembled_REs(assembled_REs):
	total_num_clusters = 0
	for chrom in assembled_REs.keys():
		total_num_clusters += len(assembled_REs[chrom])
	return total_num_clusters

def correlate_two_features(assembled_REs, feature_name_a, feature_name_b):
	list_a = []
	list_b = []
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys():
			list_a.append(assembled_REs[chrom][mycluster][feature_name_a])
			list_b.append(assembled_REs[chrom][mycluster][feature_name_b])
	return list_a, list_b

def correlate_feature_with_fc(assembled_REs, feature_name, pc = 1):
	list_a = []
	list_b = []
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys():
			list_a.append(assembled_REs[chrom][mycluster][feature_name])
			
			target_name = "index8"
			control_name = "index6"
			target_rc = assembled_REs[chrom][mycluster][target_name + "_rc"]
			control_rc = assembled_REs[chrom][mycluster][control_name + "_rc"]
			fc = (target_rc * 1.0 + pc) / (control_rc + pc)
			list_b.append(fc)
			
	return list_a, list_b
	
def convert_assembled_RE_to_single_RE_annotation(assembled_REs, summaries, target_name, control_name, target_library_size_vs_control_library_size, pc = 1):
	"""
	assembled REs: 
	{chrom: 
		{(region_start, region_end):
			{"elements":[ids]; 
			"strand": "+"; 
			target_name + "_rc":
			control_name + "_rc":
			"num_boundary_elements": value
			"5_boundary_elements":[id]; 
			"3_boundary_elements":[id];
			"I_boundary_elements":[id]; 
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	
	RE_info:
	{id:
		{
		target_name + "_rc":
		control_name + "_rc":
		target_name + "_cluster_rc":
		control_name + "_cluster_rc":
		"cluster_elements":[ids]; 
		"cluster_region":(start, end); 
		"num_boundary_elements": value;
		"expression_fc": value
		"cluster_expression_fc":max_fc
		}
	}
	
	"""
	# reorder the data so that it is indexed by id
	
	# find the summary index that the id belong to.
	RE_vs_family = find_family_for_id(summaries)
	
	RE_info = {} 
	for chrom in assembled_REs.keys():
		for mycluster in assembled_REs[chrom].keys():
			ids_in_mycluster = assembled_REs[chrom][mycluster]["elements"]
			start = mycluster[0]
			end = mycluster[1]
			expression_fc_for_ids_in_mycluster = get_expression_fc_for_RE_clusters(summaries, target_name, control_name, target_library_size_vs_control_library_size, pc, ids_in_mycluster) #{id:fc}
			for myid in ids_in_mycluster:
				
				RE_info[myid] = {}
				RE_info[myid]["cluster_region"] = (start, end)
				RE_info[myid]["cluster_elements"] = assembled_REs[chrom][mycluster]["elements"]
				RE_info[myid]["num_boundary_elements"] = assembled_REs[chrom][mycluster]["num_boundary_elements"]
				RE_info[myid]["expression_fc"] = expression_fc_for_ids_in_mycluster[myid]
				RE_info[myid]["cluster_expression_fc"] = assembled_REs[chrom][mycluster]["expression_fc_" + target_name + "_vs_" + control_name]
				RE_info[myid][target_name + "_cluster_rc"] = assembled_REs[chrom][mycluster][target_name + "_rc"]
				RE_info[myid][control_name + "_cluster_rc"] = assembled_REs[chrom][mycluster][control_name + "_rc"]
				
				summaryindex = RE_vs_family[myid]
				RE_info[myid][target_name + "_rc"] = summaries[summaryindex][myid][target_name + "_rc"]
				RE_info[myid][control_name + "_rc"] = summaries[summaryindex][myid][control_name + "_rc"]
			
	return RE_info 

def find_family_for_id(summaries):
	"""
	return {id:summary_index}
	"""
	RE_info = {}
	for i in xrange(len(summaries)):
		summary = summaries[i]
		for myid in summary.keys():
			RE_info[myid] = i
	return RE_info
	
def find_family_name(summaries):
	"""
	get the name of the family from id
	"""
	names = []
	for summary in summaries:
		myid = summary.keys()[0]
		names.append (find_name_from_id(myid))
	return names

def find_name_from_id(myid):
	p = re.compile("chr")
	my_name = p.split(myid)[0]
	return my_name
	
def find_strand_from_id(myid):
	"""
	get the strand information from id
	"""
	p = re.compile("chr")
	my_name = p.split(myid)[-1]
	plus = re.compile("\+")
	minus = re.compile("\-")
	if plus.search(my_name):
		strand = "+"
	elif minus.search(my_name):
		strand = "-"
	else:
		print my_name
		print "strand information wrong in ", myid
		sys.exit(1)
	return strand
	
	
def assign_strand_for_cluster(summaries, percentage_threshold, idset=None):
	"""
	For a cluster of ids in an assembled RE, call the collective strand
	
	summaries: [summary]
	summary:{id:{feature_name:value}}
	percentage_threshold: we pick the strang that has percentage > percentage_threshold
	
	return: "+", "-" or "I" ( for indeterminate)
	return: [strand]
	"""

	# Get the subset
	mysummaries = get_subset(summaries, idset)
	
	
	# Collect the information for strands
	strands = []
	num_positive = 0
	num_negative = 0
	total = 0
	for mysummary in mysummaries:
		for myid in mysummary.keys():
			re = mysummary[myid]["annotation"]
			strand = re.strand
			total += 1
			strands.append(strand)
			if strand == "+":
				num_positive += 1
			elif strand == "-":
				num_negative += 1
			else:
				print "Strand info is wrong for ", myid
	if num_positive * 1.0 / total >= percentage_threshold:
		strand_rep = "+"
	elif num_negative * 1.0 / total >= percentage_threshold:
		strand_rep = "-"
	else:
		strand_rep = "I"
	return strand_rep, strands
	
	
def Union_REs_in_a_family(summary, extension, strand_assignment_cutoff=0.7, idset= None):
	"""
	Motivation: Many of the incomplete sandwich structure detected by explore_sandwich_structure_for_one_RE is because this:
		RLTR4_Mm int int int int RLTR4_Mm
		the several ints in the middle are in fact pieces of a complete RLTR4_Mm_int. For unknown reason they are annotated to be different RLTR4_Mm_ints.
	Objective: union instances of REs in a single family. If we need to union REs from multiple family, we use a slightly different approach, with preprocessing the REs so that each instance is represented by (start, end, id)
	strand_assignment_cutoff: for a cluster of REs, the strand of the cluster is assigned to "I" if below the cutoff
	return:
		{chrom: {(region_start, region_end):{"elements":[ids]; "strand": "+"}}}
	"""
	
	return Union_REs_in_families([summary], extension, strand_assignment_cutoff, idset)
	
	

def Union_REs_in_families(summaries, extension, strand_assignment_cutoff=0.7, idset= None):	
	"""
	Motivation: Many of the incomplete sandwich structure detected by explore_sandwich_structure_for_one_RE is because this:
		RLTR4_Mm int int int int RLTR4_Mm
		the several ints in the middle are in fact pieces of a complete RLTR4_Mm_int or MULV_int. For unknown reason they are annotated to be different RLTR4_Mm_ints.
	Objective: union instances of REs in a single family. If we need to union REs from multiple family, we use a slightly different approach, with preprocessing the REs so that each instance is represented by (start, end, id)
	summaries: [summary]
	strand_assignment_cutoff: for a cluster of REs, the strand of the cluster is assigned to "I" if below the cutoff
	return:
		{chrom: {(region_start, region_end):{"elements":[ids]; "strand": "+"}}}
	return:
		
	"""
	# Get the subset
	mysummaries = get_subset(summaries, idset)
	
	# Separate by chrom, then resection RE families according to chrom
	my_complete_summary_by_chrom = {} #{chrom:{id:{feature_name:value}}}
	for mysummary in mysummaries:
		mysummary_by_chrom =  separate_by_chrom(mysummary) #{chrom:{id:{feature_name:value}}}
		chroms = mysummary_by_chrom.keys()
		for chrom in chroms:
			if chrom not in my_complete_summary_by_chrom.keys():
				my_complete_summary_by_chrom[chrom] = {}
			my_complete_summary_by_chrom[chrom].update(mysummary_by_chrom[chrom])
			#print chrom, len(my_complete_summary_by_chrom[chrom])
	
	
	# Cluster the REs: {chrom: {(region_start, region_end):{"elements":[ids]; "strand": "+"; }}}
	clustered_REs = {}
	chroms = my_complete_summary_by_chrom.keys()
	for chrom in chroms:
		REs_by_chrom = my_complete_summary_by_chrom[chrom] #{id:{feature_name:value}}
		REs_data = []
		for myid in REs_by_chrom.keys():
			myre = REs_by_chrom[myid]
			re_element = myre['annotation']
			start, end = re_element.genoStart, re_element.genoEnd
			REs_data.append((start, end, myid)) # [((start, end, myid))]
		assembles = Utility_extended.union_with_trace(REs_data, extension) # {(start, end):[(start, end, id)]}
		clustered_REs[chrom] = {}
		for region in assembles.keys():
			clustered_REs[chrom][region] = {}
			ids = [item[2] for item in assembles[region]]# item[2] is id
			clustered_REs[chrom][region]["elements"] =  ids
			# strand_rep = "+", "-", "I"
			strand_rep, strands	 = assign_strand_for_cluster(mysummaries, strand_assignment_cutoff, ids)
			clustered_REs[chrom][region]["strand"] = strand_rep	
	return clustered_REs	
		
	
def separate_by_chrom(re_elements):
	"""
	re_elements: {id:{feature_name:value}}
	return {chrom:{id:{feature_name:value}}}
	"""
	re_elements_by_chrom = {}
	for myid in re_elements.keys():
		element = re_elements[myid]['annotation']
		if element.chrom not in re_elements_by_chrom.keys():
			re_elements_by_chrom[element.chrom] = {}
		re_elements_by_chrom[element.chrom][myid] = re_elements[myid]
	return re_elements_by_chrom

def select_by_age(elements, age_threshold):
	"""
	select elements that are younger than certain age_threshold
	elements: {id:{feature_name:value}}
	return: {id:{feature_name:value}}
	"""
	selected_elements = {}
	for myid in elements.keys():
		element = elements[myid]['annotation']
		if element.age <= age_threshold:
			selected_elements[myid] = elements[myid]
	return selected_elements
	

	
def get_shared(de_members1, de_members2):
	"""
	de_members: {reClass:{reFamily:reName:[id]}}}
	shared_members:  {reClass:{reFamily:reName:[id]}}}
	"""
	shared_members ={}
	myClasses = list(set(de_members1.keys()) & set(de_members2.keys()) )
	for reClass in myClasses:
		shared_members[reClass] = {}
		myFamilies = list(set(de_members1[reClass].keys()) & set(de_members2[reClass].keys()) )
		for reFamily in myFamilies:
			shared_members[reClass][reFamily] = {}
			myNames = list(set(de_members1[reClass][reFamily].keys()) & set(de_members2[reClass][reFamily].keys()) )
			for reName in myNames:
				shared_members[reClass][reFamily][reName] = list(set(de_members1[reClass][reFamily][reName]) & set(de_members2[reClass][reFamily][reName]) )
	return shared_members

def plot_bar_charts(dic, outfilename):
	"""
	dic = {name: value}
	"""
	N = len(dic)
	x = np.arange(1, N+1)
	sorted_tuple_list = sorted(dic.items(), key=itemgetter(1), reverse=True)
	y = [ item[1] for item in  sorted_tuple_list]
	labels = [item[0] for item in sorted_tuple_list ]
	width = 1
	plt.clf()
	bar1 = plt.bar( x, y, width, color="y" )
	plt.ylabel( '' )
	
	plt.tick_params(axis='x', which='major', labelsize=6)
	plt.tick_params(axis='x', which='minor', labelsize=6)
	plt.xticks(x + width/2.0, labels, rotation=45)
	plt.savefig(outfilename, format="eps")
	plt.close()
    
def plot_pie_charts(dic, outfilename):
	"""
	dic = {name: value}
	"""
	sorted_tuple_list = sorted(dic.items(), key=itemgetter(1), reverse=True)
	y = [ item[1] for item in  sorted_tuple_list]
	labels = [item[0] for item in sorted_tuple_list ]
	plt.figure()
	plt.suptitle(outfilename)
	wedges, labels = plt.pie(y, labels=labels, explode=None, shadow=False)
	plt.axis('equal')
	plt.savefig(outfilename, format="eps")
	plt.close()

def output_simple_dic(dic, outfilename):
	out = open(outfilename, "w")
	sorted_tuple_list = sorted(dic.items(), key=itemgetter(1), reverse=True)
	for item in sorted_tuple_list:
		myid = item[0]
		value = item[1]
		out.write(myid + "\t" + str(value) + "\n")
	out.close()
    
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
	
	##load the annotation files and add annotation to features if that is not done
	#print "\n\nExamine features and load annotation if necessary"
	#index = 0
	#for reClass in re_tree.keys():
		#for reFamily in re_tree[reClass].keys():
			#for reName in re_tree[reClass][reFamily]:
				#summary_file_name = opt.summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				#assert( Utility_extended.fileExists(summary_file_name) == 1)
				#inf = open(summary_file_name, 'rb')
				#summary = pickle.load(inf)
				#inf.close()
				#index += 1
				#feature_names = get_feature_names(summary) 
				#if index <2: #only print once
					#print "The features stored in the ", summary_file_name, " are ",  feature_names
				#if 'annotation' not in feature_names:
					#print "undating %s" %summary_file_name
					#add_annotation_to_features(opt.RE_pkl_file_location, summary, opt.summary_name, reClass, reFamily, reName)
					
	fc = 4
	min_rpkm = 0.5
	min_rc = 5
	
	target_name = "index8"
	control_name = "index6"
	print "\n\n", target_name, " vs ",  control_name
	
	target_library_size_vs_control_library_size = get_target_library_size_vs_control_library_size(re_tree, opt.summary_name, target_name, control_name)
	
	if target_library_size_vs_control_library_size == -1:
		print "There is not a single locus with reads in both %s and %s" %(target_name, control_name)
		exit(1)
	
	#(de_members_86, de_enrichment, total_instances) = get_all_upregulated_REs(re_tree, opt.summary_name, target_name, control_name, target_library_size_vs_control_library_size, fc, min_rpkm, min_rc, pc = 1)
	
	#out_file_name = target_name + "_vs_" + control_name + "_fc" + str(fc) + "DE_enrichment_ranking.txt"
	#mylist = rank_enrichment(de_members_86, de_enrichment, total_instances, out_file_name)
	
	#for i in xrange(10):
		#reClass = mylist[i][0]
		#reFamily = mylist[i][1]
		#reName = mylist[i][2]
		
		#summary_file_name = opt.summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
		#inf = open(summary_file_name, 'rb')
		#summary = pickle.load(inf)
		#if i == 0:
			#print get_feature_names(summary)
		#inf.close()
		
		## output instances of differentially expressed regions
		#out_file_name = target_name + "_vs_" + control_name + "_DE_elements_fc" + str(fc) + "_on_" + "_".join([reClass, reFamily, reName]) + ".txt"
		#output_instances(de_members_86, summary, reClass, reFamily, reName, out_file_name)
	
		
		#name = "_".join([reClass, reFamily, reName])
		#get_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
		#get_log_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
		#get_expre_vs_fc(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
	
	
	#output = open("86_derepressed_instances.pkl", 'wb')
	#pickle.dump(de_members_86, output)
	#output.close()
	
	#output = open("RE_instance_counts.pkl", 'wb')
	#pickle.dump(total_instances, output)
	#output.close()
	
	print "\n\nLoading 86_derepressed_instances.pkl"
	inf = open("86_derepressed_instances.pkl", 'rb')
	de_members_86 = pickle.load(inf)
	inf.close()
	print "Loading RE_instance_counts.pkl"
	inf = open("RE_instance_counts.pkl", 'rb')
	total_instances = pickle.load(inf)
	inf.close()
	
	class_level_enrichment_counts = get_class_level_enrichment_counts(de_members_86)
	output_simple_dic(class_level_enrichment_counts, "86_class_level_enrichment_counts.txt")
	non_zero_class_level_enrichment_counts = get_entries_above_cutoff(class_level_enrichment_counts, cutoff=0)
	output_simple_dic(non_zero_class_level_enrichment_counts, "86_non_zero_class_level_enrichment_counts.txt")
	print sorted(non_zero_class_level_enrichment_counts.items(), key=itemgetter(1), reverse=True)
	plot_bar_charts(class_level_enrichment_counts, "86_class_level_enrichment_counts_barplot.eps")
	#plot_pie_charts(class_level_enrichment_counts, "class_level_enrichment_counts_pieplot.eps")
	#plot_bar_charts(class_level_enrichment_counts, "86_class_level_enrichment_counts_barplot.eps")
	plot_pie_charts(non_zero_class_level_enrichment_counts, "86_non_zero_class_level_enrichment_counts_pieplot.eps")
	
	class_level_enrichment_ratios, above_cutoff_class_level_enrichment_ratios = get_class_level_enrichment_ratios(de_members_86, total_instances, cutoff=10)
	print above_cutoff_class_level_enrichment_ratios
	output_simple_dic(class_level_enrichment_ratios, "86_class_level_enrichment_ratios.txt")
	output_simple_dic(above_cutoff_class_level_enrichment_ratios, "86_above_cutoff_class_level_enrichment_ratios.txt")
	plot_bar_charts(above_cutoff_class_level_enrichment_ratios, "86_above_cutoff_class_level_enrichment_ratios_barplot.eps")
	#plot_pie_charts(above_cutoff_class_level_enrichment_ratios, "86_above_cutoff_class_level_enrichment_ratios_pieplot.eps")
	
	#print de_members_86.keys()
	reClass = 'LTR'
	#reClass = 'SINE'
	family_level_enrichment_counts = get_family_level_enrichment_counts(de_members_86, reClass)
	output_simple_dic(family_level_enrichment_counts, "86_" + reClass+ "_family_level_enrichment_counts.txt")
	non_zero_family_level_enrichment_counts = get_entries_above_cutoff(family_level_enrichment_counts, cutoff=0)
	output_simple_dic(non_zero_family_level_enrichment_counts, "86_" + reClass + "_non_zero_family_level_enrichment_counts.txt")
	print sorted(non_zero_family_level_enrichment_counts.items(), key=itemgetter(1), reverse=True)
	plot_bar_charts(family_level_enrichment_counts, "86_" + reClass+ "_family_level_enrichment_counts_barplot.eps")
	plot_pie_charts(non_zero_family_level_enrichment_counts, "86_" + reClass + "_non_zero_family_level_enrichment_counts_pieplot.eps")
	
	family_level_enrichment_ratios, above_cutoff_family_level_enrichment_ratios = get_family_level_enrichment_ratios(de_members_86, reClass, total_instances, cutoff = 0)
	print sorted(above_cutoff_family_level_enrichment_ratios.items(), key=itemgetter(1), reverse=True)
	output_simple_dic(family_level_enrichment_ratios, "86_" + reClass + "_family_level_enrichment_ratios.txt")
	plot_bar_charts(family_level_enrichment_ratios, "86_" + reClass + "_family_level_enrichment_ratios_barplot.eps")

	
	
	
	reClass = 'LTR'
	reFamily = "ERV1"
	name_level_enrichment_counts = get_name_level_enrichment_counts(de_members_86, reClass, reFamily)
	output_simple_dic(name_level_enrichment_counts, "86_" + reClass+ "_" + reFamily + "_name_level_enrichment_counts.txt")
	non_zero_name_level_enrichment_counts = get_entries_above_cutoff(name_level_enrichment_counts, cutoff=0)
	print sorted(non_zero_name_level_enrichment_counts.items(), key=itemgetter(1), reverse=True)
	output_simple_dic(non_zero_name_level_enrichment_counts, "86_" + reClass+ "_" + reFamily + "_non_zero_name_level_enrichment_counts.txt")
	plot_bar_charts(name_level_enrichment_counts, "86_" + reClass+ "_" + reFamily + "_name_level_enrichment_counts_barchart.eps")
	name_level_enrichment_ratios = get_name_level_enrichment_ratios(de_members_86, reClass, reFamily, total_instances)
	output_simple_dic(name_level_enrichment_ratios, "86_" + reClass+ "_" + reFamily + "_name_level_enrichment_ratios.txt")
	non_zero_name_level_enrichment_ratios = get_entries_above_cutoff(name_level_enrichment_ratios, cutoff=0)
	output_simple_dic(non_zero_name_level_enrichment_ratios, "86_" + reClass + "_" + reFamily + "_non_zero_name_level_enrichment_ratios.txt")
	print sorted(non_zero_name_level_enrichment_ratios.items(), key=itemgetter(1), reverse=True)
	plot_bar_charts(non_zero_name_level_enrichment_counts, "86_" + reClass + "_" + reFamily + "_non_zero_name_level_enrichment_counts_barchart.eps")
	plot_pie_charts(non_zero_name_level_enrichment_counts, "86_" + reClass + "_" + reFamily + "_non_zero_name_level_enrichment_counts_piechart.eps")
	plot_bar_charts(non_zero_name_level_enrichment_ratios, "86_" + reClass + "_" + reFamily + "_non_zero_name_level_enrichment_ratios_barchart.eps")
	
	
	
	#target_name = "index9"
	#control_name = "index6"
	#print target_name, " vs ",  control_name
	#target_library_size_vs_control_library_size = get_target_library_size_vs_control_library_size(re_tree, opt.summary_name, target_name, control_name)
	
	#if target_library_size_vs_control_library_size == -1:
	#	print "There is not a single locus with reads in both %s and %s" %(target_name, control_name)
	#	exit(1)
		
		
	##de_members = {} #{reClass:{reFamily:reName:[id]}}}
	##de_enrichment = {}  #{reClass:{reFamily:reName:enrichment}}}
	##total_instances = {} #{reClass:{reFamily:reName:count}}}
	#(de_members_96, de_enrichment, total_instances) = get_all_upregulated_REs(re_tree, opt.summary_name, target_name, control_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc=15, pc = 1)
	#out_file_name = target_name + "_vs_" + control_name + "_fc" + str(fc) + "DE_enrichment_ranking.txt"
	#mylist = rank_enrichment(de_members_96, de_enrichment, total_instances, out_file_name)
	
	#for i in xrange(10):
		#reClass = mylist[i][0]
		#reFamily = mylist[i][1]
		#reName = mylist[i][2]
		
		#summary_file_name = opt.summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
		#inf = open(summary_file_name, 'rb')
		#summary = pickle.load(inf)
		#inf.close()
		
		## output instances of differentially expressed regions
		#out_file_name = target_name + "_vs_" + control_name + "_DE_elements_fc" + str(fc) + "_on_" + "_".join([reClass, reFamily, reName]) + ".txt"
		#output_instances(de_members_96, summary, reClass, reFamily, reName, out_file_name)
	
		
		##name = "_".join([reClass, reFamily, reName])
		##get_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
		##get_log_fc_histogram(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
		##get_expre_vs_fc(summary, target_name, control_name, name, target_library_size_vs_control_library_size, pc = 1)
	
	#shared_instances = get_shared(de_members_86, de_members_96)
	#output_all_instances(shared_instances, opt.summary_name, "86_vs_96_fc4_shared.txt")
	
	#target_name = "index12"
	#control_name = "index9"
	#print target_name, " vs ",  control_name
	#target_library_size_vs_control_library_size = get_target_library_size_vs_control_library_size(re_tree, opt.summary_name, target_name, control_name)
	
	#if target_library_size_vs_control_library_size == -1:
	#	print "There is not a single locus with reads in both %s and %s" %(target_name, control_name)
	#	exit(1)
	
	#(de_members, de_enrichment, total_instances) = get_all_upregulated_REs(re_tree, opt.summary_name, target_name, control_name, fc, target_library_size_vs_control_library_size, min_rpkm, min_rc=15, pc = 1)
	#out_file_name = target_name + "_vs_" + control_name + "_fc" + str(fc) + "DE_enrichment_ranking.txt"
	#mylist = rank_enrichment(de_members, de_enrichment, total_instances, out_file_name)
	
	#for i in xrange(10):
		#reClass = mylist[i][0]
		#reFamily = mylist[i][1]
		#reName = mylist[i][2]
		
		#summary_file_name = opt.summary_name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
		#inf = open(summary_file_name, 'rb')
		#summary = pickle.load(inf)
		#inf.close()
		
		## output instances of differentially expressed regions
		#out_file_name = target_name + "_vs_" + control_name + "_DE_elements_fc" + str(fc) + "_on_" + "_".join([reClass, reFamily, reName]) + ".txt"
		#output_instances(de_members, summary, reClass, reFamily, reName, out_file_name)
		
		
		
		
if __name__ == "__main__":
	main(sys.argv)
