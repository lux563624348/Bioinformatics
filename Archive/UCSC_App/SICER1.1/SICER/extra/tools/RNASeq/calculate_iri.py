#!/usr/bin/env python
"""
Function: Calculate and characterize the IRI of one state (eg, resting)

The IRI calculated is slightly different from AnalyzeRNASeqPASeqChIPSeq.py:
IRI: Previously I used (rpkm+pc/rpkm+pc)
    Now I use ((rc +pc )/length) /((rc + pc)/length) 
    
required as input: 

one summary pickle data structures, obtained from get_non_strandspecific_read_count_on_ExonsIntrons.py/get_strandspecific_read_count_on_ExonsIntrons.py

each: {entrez_id:{feature_name:value}}

	summary: {entrezID:{attribute:value}}
				summary[entrez_id] = {}
				(summary[entrez_id])["merged_exons_rc"] = merged_exons_rc
				(summary[entrez_id])["merged_exon_RPKM"] = merged_exon_RPKM
				(summary[entrez_id])["merged_exons_total_length"] = merged_exons_total_length
				(summary[entrez_id])["shared_exons_rc"] = shared_exons_rc
				(summary[entrez_id])["shared_exon_RPKM"] = shared_exon_RPKM
				(summary[entrez_id])["shared_exons_total_length"] = shared_exons_total_length
				(summary[entrez_id])["shared_introns_rc"] = shared_introns_rc
				(summary[entrez_id])["shared_intron_RPKM"] = shared_intron_RPKM
				(summary[entrez_id])["shared_introns_total_length"] = shared_introns_total_length
				(summary[entrez_id])["merged_transcript_rc"] = merged_transcript_rc
				(summary[entrez_id])["merged_transcript_RPKM"] = merged_transcript_RPKM
				(summary[entrez_id])["merged_transcript_length"] = merged_transcript_length




"""

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
import matplotlib
import copy

try:
   import cPickle as pickle
except:
   import pickle

import rpy2.robjects as robjects
r = robjects.r
from rpy2.robjects.packages import importr

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/JunZhu/modules')

import AnalyzeRNASeqPASeqChIPSeq as master
import AnalyzeRNASeqChIPSeq
import GenomeData
import Utility_extended

#font = {'family' : 'normal',
        #'weight' : 'bold',
        #'size'   : 14}
#matplotlib.rc('font', **font)

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")


def calculate_expression_iri(expression_summary_dic, pc = 1):
	"""
	If intron length = 0 or the gene is not expressed at all: iri = 0
	pc is the read count pseudo-count
	Else:
		iri[myid] = shared_intron_rpkm/shared_exon_rpkm 
	retains everything and add the "iri" to the dic
	
	return: {entrez_id:{feature_name:value}}
	"""
	iri ={}
	
	for myid in expression_summary_dic.keys():
		total_introns_length = float(expression_summary_dic[myid]["shared_introns_total_length"])
		total_exons_length = expression_summary_dic[myid]["shared_exons_total_length"]
		if total_introns_length == 0:
			iri[myid] = 0.0
		elif total_exons_length == 0:
			assert expression_summary_dic[myid]["shared_exons_rc"] < 1
		else:
			shared_intron_rpkm = (float(expression_summary_dic[myid]["shared_introns_rc"]) + pc)/float(expression_summary_dic[myid]["shared_introns_total_length"])
			shared_exon_rpkm = (float(expression_summary_dic[myid]["shared_exons_rc"]) + pc)/float(expression_summary_dic[myid]["shared_exons_total_length"])
			iri[myid] = shared_intron_rpkm/shared_exon_rpkm
			# If the gene is not expressed at all, the iri = 0
			if expression_summary_dic[myid]["shared_introns_rc"] == 0 and expression_summary_dic[myid]["shared_exons_rc"] == 0:
				iri[myid] = 0.0 
	
	new_expression_summary_dic = master.add_feature(expression_summary_dic, "iri", iri)
	return new_expression_summary_dic
	
def removed_misannotated_ids_from_IRI(summary, max_iri = 1):
	"""
	Remove ids whose intron length is 0,  and whose iri > 1
	
	Returns the ids
	"""
	cluster0 = []
	for myid in summary.keys():
		if summary[myid]["shared_introns_total_length"] > 0:
			if summary[myid]["iri"]<=max_iri:
				cluster0.append(myid)
	return cluster0

def get_expressed_genes(summary,  minRPKM):
	myids=[]
	for myid in summary.keys():
		RPKM = (summary[myid])["shared_exon_RPKM"]
		if RPKM >= minRPKM:
			myids.append(myid)
	if myid == 1837:
		print
		print (summary[myid])["shared_exon_RPKM"]
		print 
	return myids	
	
def characterize_IRI(summary, name, pc=0.00001):
	"""
	Draw IRI histograms
	pc is for log transformation
	
	"""
	#font = {'family':'normal', 'weight':'bold', 'size':18}
	#matplotlib.rc('font', **font)
	keys = sorted(summary.keys())
	iri_list = [summary[id]["iri"] for id in keys]
	print "IRI before log2: mean", numpy.average(iri_list), "    median:",  numpy.median(iri_list)
	(num_iri_0, num_iri_greater_than_1) = master.iri_component(iri_list)
	print "There are (%f, %f) genes with iri (=0, >=1)" %(float(num_iri_0)/len(iri_list), float(num_iri_greater_than_1)/len(iri_list)) 
	log_iri_list = numpy.array([log(i + pc, 2) for i in iri_list])
	print "IRI after log2:  mean ",  numpy.average(log_iri_list), "    median:", numpy.median(log_iri_list)
	
	#Histogram
	plt.clf()
	plt.hist(log_iri_list, bins=100, color='darkgreen', normed=True, alpha=0.5,label="")
	myfontsize=28
	plt.title("IRI histogram", fontsize=myfontsize)
	plt.xlabel("log(IRI)", fontsize=myfontsize)
	plt.ylabel("Frequency", fontsize=myfontsize)
	#plt.legend()
	plt.legend(loc = 'upper left')
	if name !="":
		plt.savefig( name + "_iri_hist.eps", format="eps")
		plt.savefig( name + "_iri_hist.png", format="png")
	plt.close()
	
	#Boxplot
	myboxplot = plt.boxplot(log_iri_list, sym='')
	plt.setp(myboxplot['medians'], color = 'black')
	plt.setp(myboxplot['boxes'], color = 'black' )
	plt.setp(myboxplot['whiskers'], color = 'black')
	plt.setp(myboxplot['fliers'], color = 'black')
	if name !="":
		title = name + "_iri_boxplot"
		plt.savefig(title + ".eps", format="eps")	
	
	return iri_list
		
def rank_by_iri(summary):
	"""
	iri high is at the top
	Return: [id]
	"""
	iri_id_list_w_iri = [(myid, summary[myid]["iri"]) for myid in summary.keys()] #[(id, iri)]
	iri_id_list_w_iri.sort(key=itemgetter(1), reverse = True)	
	iri_id_list = [item[0] for item in iri_id_list_w_iri]
	return iri_id_list
	
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--all_feature_pickle_file", action="store", type="string", dest="summary", help="summary of  RNASeq features in a pickle file", metavar="<file>")
	parser.add_option("-n", "--name", action="store", type="string", dest="name", help="name for the data set", metavar="<file>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	# load RNA features
	if 	Utility_extended.fileExists(opt.summary): #directly read the pickle file if it is available
		expression_summary = pickle.load(open(opt.summary, 'rb'))
		print "There are %d entrez IDs in %s" %(len(expression_summary.keys()), opt.summary)
		print "From the pickle there are %d features" % master.get_num_features(expression_summary)
		print master.get_feature_names(expression_summary)
	else:	
		print opt.summary, " is not found";
		sys.exit(1)
		
	"""
	Calculate IRI
	"""
	iri_pc = 1 # read count pc
	# expressed_summary is a dic
	expression_summary = calculate_expression_iri(expression_summary, iri_pc)
	print "After adding iri. There are %d features in %s" % (master.get_num_features(expression_summary), opt.summary)
		
	"""
	characterize IRI
	"""
	print "\nRemove genes with IRI>=1 and genes without intron"
	
	# Remove genes with likely misannotations, remove genes with 0 introns
	max_iri = 1
	ids_misannotation_removed = removed_misannotated_ids_from_IRI(expression_summary, max_iri)
	print "\nAfter removing mis-annotations and genes with no introns, there are %d IDs left for downstream analysis" %len(ids_misannotation_removed)
	
	# Find expressed genes
	# cutoff 1 is better than cutoff 0.5 for histone mark iri results, cutoff 3 is better in H3K36me3 iri top 1000 vs bottom 1000
	cutoff = 1
	print "\nRemove genes with resting expression rpkm <1"
	expressed_ids = get_expressed_genes(expression_summary, cutoff)
	print "\nWith requirement of resting RPKM above ", cutoff, ", ",  len(expressed_ids), " entrez IDs are left for downstream analysis"

	#Combine the two requirement
	print "\nRemove genes with resting expression rpkm <1 and with misannotation"
	subset_ids = list (set(ids_misannotation_removed) & set(expressed_ids))
	print "\nAfter removing mis-annotations,genes with no introns and genes not expressed in %s, there are %d IDs left for downstream analysis" %(opt.summary, len(subset_ids))
	
	expressed_summary = master.subset(subset_ids, expression_summary)
	
	print "\nCharacterize IRI"
	iri_pc = 0.001
	characterize_IRI(expressed_summary, opt.name, iri_pc)	
	
	#Rank the genes according to iri, and select the top 1000 and bottom 1000 for analysis
	number_of_bins = 6
	print "\n\nTop 1000 vs bottom 1000" 
	ids_ranked_by_iri = rank_by_iri(expressed_summary)
	number_of_genes = 1000
	iri_top_ids = ids_ranked_by_iri[:number_of_genes]
	iri_bottom_ids = ids_ranked_by_iri[-1*number_of_genes:]
	#subset
	iri_top1000 = master.subset(iri_top_ids, expressed_summary)
	iri_bottom1000 = master.subset(iri_bottom_ids, expressed_summary)
	#bin the ids [[id]]
	iri_top_ids_binned = master.bin_target_genes_by_feature(iri_top_ids, expressed_summary, 'shared_exon_RPKM', number_of_bins)
	iri_bottom_ids_binned = master.bin_target_genes_by_feature(iri_bottom_ids, expressed_summary, 'shared_exon_RPKM', number_of_bins)
	sets = (iri_top_ids_binned, iri_bottom_ids_binned)
	sets_name = opt.name + "_iri_top" + str(number_of_genes) + "_bottom" + str(number_of_genes)
	
	pc = 0.001
	# Validation of partitioning on iri
	print "\nLog iri"
	feature = dict((k, log(expressed_summary[k]["iri"]+pc, 2)) for k in expressed_summary.keys())
	title = "iri_" + sets_name + "_boxplot"
	master.compare_feature_on_sets_of_binned_genes(sets, feature, title)
	
	# Validation of partitioning on expression
	print "\nLog expression"
	feature = dict((k, log(expressed_summary[k][ "shared_exon_RPKM"], 2)) for k in expressed_summary.keys())
	title = "expression_" + sets_name + "_boxplot"
	master.compare_feature_on_sets_of_binned_genes(sets, feature, title)
	
	print "\n\nCompare the exon length and transcript length between iri top1000 and bottom 1000 genes"
	# The latest finding is that top1000 has much shorter genes that bottom1000, which is mainly due to the shorter introns. The concern is that the difference in IRI could be due to coverage, ie, the long introns will exhibit lower intronic read density due to lack of coverage. To address this concern, we can require that the intronic RPKM to be bigger than a threshold. The result confirmed the finding.
	myids = []
	intron_RPKM_threshold = 0.1
	for myid in expressed_summary.keys():
		if expressed_summary[myid]["shared_intron_RPKM"] >= intron_RPKM_threshold:
			myids.append(myid)
	intron_expressed_summary= master.subset(myids, expressed_summary)
	print "Compare the exon length and transcript length between iri top1000 and bottom 1000 genes which are expressed in intron"
	for feature_name in ["merged_exons_total_length", 'merged_transcript_length']:
		title =  feature_name + "_intron_expressed_" + sets_name + "_boxplot"
		feature_dic= dict((k, log(expressed_summary[k][feature_name],2)) for k in intron_expressed_summary.keys())
		master.compare_feature_on_sets_of_binned_genes(sets, feature_dic, title)
	
if __name__ == "__main__":
	main(sys.argv)	
	
