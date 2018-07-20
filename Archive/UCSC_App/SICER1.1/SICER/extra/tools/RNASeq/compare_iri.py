#!/usr/bin/env python
"""
Compare the IRI of multiple states (eg, resting vs activated)

The iri_change vs expression_change: the anti-correlation can be perceived as by definition. 
Solution: use splice-junction reads to gauge expression change. The splice-junction reads are not used in calculating IRI.


required as input: 

multiple summary pickle data structures, each of which is obtained from get_non_strandspecific_read_count_on_ExonsIntrons.py/get_strandspecific_read_count_on_ExonsIntrons.py

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
import calculate_iri
import GenomeData
import Utility_extended

#font = {'family' : 'normal',
        #'weight' : 'bold',
        #'size'   : 14}
#matplotlib.rc('font', **font)

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")


def get_names(summary_file):
	"""
	Read the file that contains all the name and file_name of the libraries
	
	The format of the file is expected to be:
	name	file_name
	name	file_name
	
	The order of the name is useful information, the order information is stored in names
	
	Return:
	[name]
	{name:file_name}
	"""
	
	of = open(summary_file, "r")
	names = []
	summary_names = {}
	for line in of:
		if not re.match(comment, line):
			line = line.strip()
			sline = line.split()
			if len(sline) > 0:
				name = sline[0]
				file_name = sline[1]
				names.append(name)
				summary_names[name] = file_name
	of.close()
	return names, summary_names
	
def removed_misannotated_ids_from_IRI(my_summaries, max_iri = 1):
	"""
	
	my_summaries: {name:summary}
	
	For my_summaries with 2 entries:
	
	scatter plot of resting vs activated produce 4 clusters. Three of the 4 have outrageous IRI in either of the cell type and the IRI is strongly dependent on pseudo count involved in defining IRI (intron RPKM + pc)/(exon RPKM + pc) . suggesting 1) mis-annotation and furthermore 2) for those genes, the exon read counts are 0, whereas the intron read counts are non-zero. For all meaningful analysis along the IRI line, we need to get rid of those genes. Perhaps it is useful to revisit the two differential IRI outlier clusters at a later stage. This module is for use of finding the IRI main cluster.
	
	The delination of the clusters are based on the IRI scatter plot:
	cluster 0 is the main cluster, independent of the pseudo count.
	cluster 1,2 IDs to be determined
	
	Also remove ids whose intron length is 0
	
	Returns a list of legitimate IDs.
	"""
	
	first_summary = my_summaries[(my_summaries.keys()[0])]
	all_ids = first_summary.keys()
	
	illegitimate_ids = []
	for my_summary in my_summaries.values():
		for myid in my_summary.keys():
			if my_summary[myid]["shared_introns_total_length"] < 1 or my_summary[myid]["iri"] > max_iri:
				illegitimate_ids.append(myid)
	illegitimate_ids = list(set(illegitimate_ids))
	legitimate_ids = list(set(all_ids) - set(illegitimate_ids))		
	return legitimate_ids


def get_expressed_genes(summaries,  minRPKM, flag):
	"""
	summaries:{name:summary}
	flag:[name]
	
	Three styles, 
	1) require that the genes are expressed in all states flag = my_summary.keys() 
	2) In a particular set of states, flag = subset of my_summary.keys() 
	3) In any one state,  flag =[-1] 
	
	Returns [id]
	"""
	first_summary = summaries[(summaries.keys()[0])]
	all_ids = first_summary.keys()
	
	myids=[]
	if flag[0] == -1: #expressed in any one state
		for myid in all_ids:
			for name in summaries.keys():
				mysummary = summaries[name]
				if myid not in mysummary.keys():
					break
				else:
					if mysummary[myid]["shared_exon_RPKM"] >= minRPKM:
						myids.append(myid)
						break
					
	else:#expressed in designated states
		for myid in all_ids:
			yes = 1
			for name in flag:
				mysummary = summaries[name]
				if myid not in mysummary.keys() or mysummary[myid]["shared_exon_RPKM"] < minRPKM:
					yes = 0
					break
			if yes == 1:
				myids.append(myid)
	return myids	

def get_subset(summaries, subset_ids):
	new_summaries = {}
	for name in summaries.keys():
		summary = summaries[name]
		new_summary = {}
		for myid in subset_ids:
			if myid in summary.keys():
				new_summary[myid] = copy.deepcopy(summary[myid])
		new_summaries[name] = new_summary
	return new_summaries
	
def characterize_IRI(summaries, names, title, pc = 0.00001, subset_ids = None):
	"""
	summaries: {name:summary}
	names:[name] has order information
	characterize IRI by boxplots
	
	if the number of libraries is 2: histograms
	
	pc is for log transformation in Characterize_IRI
	
	"""
	if subset_ids is not None:
		newsummaries = get_subset(summaries, subset_ids)
	else:
		newsummaries = summaries
	
	log_iris = {} #{name:log_iri}
	for name in newsummaries.keys():
		newsummary = newsummaries[name]
		iri_list = [newsummary[id]["iri"] for id in newsummary.keys()]
		log_iri = numpy.array([log(i + pc, 2) for i in iri_list])
		log_iris[name] = log_iri
	

	#Boxplot
	log_iri_list = []
	for name in names:
		log_iri_list.append(log_iris[name])
		
	plt.clf()
	myboxplot = plt.boxplot(log_iri_list, sym='')
	plt.setp(myboxplot['medians'], color = 'black')
	plt.setp(myboxplot['boxes'], color = 'black')
	plt.setp(myboxplot['whiskers'], color = 'black')
	plt.setp(myboxplot['fliers'], color = 'black')
	plt.savefig(title + "_boxplot.eps", format="eps")	
	plt.close()
	
	colors = ['darkgreen', 'darkorange', 'darkpink', 'darkblue']
	if len(names) == 2:
		#Histogram
		plt.clf()
		index = 0
		for name in names:
			plt.hist(log_iris[name], bins=100, color=colors[index], normed=True, alpha=0.5,label=name)
			index += 1
		myfontsize=28
		plt.title("IRI histogram", fontsize=myfontsize)
		plt.xlabel("log(IRI)", fontsize=myfontsize)
		plt.ylabel("Frequency", fontsize=myfontsize)
		plt.legend(loc = 'upper left')
		plt.savefig( title + "_iri_hist.eps", format="eps")
		plt.savefig( title + "_iri_hist.png", format="png")
		plt.close()
		print "	P-value associated with ranking of the two sets of observations using MWU test: ", scipy.stats.mannwhitneyu( log_iri_list[0], log_iri_list[1])[1]	
	

def get_subset_genes_by_iri_changes(summary_resting, summary_activated, change_type, fc_threshold, minRPKM, iri_pc=0.001):
	"""
	change_type: "up", "down", "no-change"
	iri_pc is for iri fc
	minRPKM: as long as one side meets the requirement
	fc_threshold is activated_iri/resting_iri
	Only retain genes with shared-exon expression
	
	Return: [ids]
	"""
	common_ids = list(set(summary_resting.keys()) & set(summary_activated.keys()))
	
	myids = []
	for myid in common_ids:
		resting_RPKM = float(summary_resting[myid]["shared_exon_RPKM"]) #shared exon
		activated_RPKM = float(summary_activated[myid]["shared_exon_RPKM"]) #shared exon
		if resting_RPKM > minRPKM or activated_RPKM > minRPKM:
			iri_fc = (summary_activated[myid]["iri"] + iri_pc) / (summary_resting[myid]["iri"] + iri_pc)
			if change_type == "down":
				assert (fc_threshold<1)
				if iri_fc <= fc_threshold:
					myids.append(myid)
			elif change_type == "up":
				assert ( fc_threshold > 1)
				if iri_fc >= fc_threshold:
					myids.append(myid)
			elif change_type == "no-change":
				assert ( fc_threshold > 1)
				if iri_fc < fc_threshold and iri_fc > 1.0/fc_threshold:
					myids.append(myid)
	return myids

def get_subset_genes_by_expression_changes(summary_resting, summary_activated, change_type, fc_threshold, minRPKM, pc=0.001):
	"""
	change_type: "up", "down", "no-change"
	pc is for expression_fc
	minRPKM: as long as one side meets the requirement
	fc_threshold is activated/resting
	Only retain genes with shared-exon expression
	
	Return: [ids]
	"""
	common_ids = list(set(summary_resting.keys()) & set(summary_activated.keys()))
	
	myids = []
	for myid in common_ids:
		resting_RPKM = float(summary_resting[myid]["shared_exon_RPKM"]) #shared exon
		activated_RPKM = float(summary_activated[myid]["shared_exon_RPKM"]) #shared exon
		if resting_RPKM > minRPKM or activated_RPKM > minRPKM:
			expression_fc = (activated_RPKM + pc) / (resting_RPKM + pc)
			if change_type == "down":
				assert (fc_threshold<1)
				if expression_fc <= fc_threshold:
					myids.append(myid)
			elif change_type == "up":
				assert ( fc_threshold > 1)
				if expression_fc >= fc_threshold:
					myids.append(myid)
			elif change_type == "no-change":
				assert ( fc_threshold > 1)
				if expression_fc < fc_threshold and expression_fc > 1.0/fc_threshold:
					myids.append(myid)
	return myids

def get_expression_changes (summary_resting, summary_activated, id_set, pc=0.001):
	"""
	returns [(id, expression_fc)]
	"""
	mylist = []
	for myid in id_set:
		resting_RPKM = float(summary_resting[myid]["shared_exon_RPKM"]) #shared exon
		activated_RPKM = float(summary_activated[myid]["shared_exon_RPKM"]) #shared exon
		expression_fc = (activated_RPKM + pc) / (resting_RPKM + pc)
		mylist.append([myid, expression_fc])
	return mylist
	
	
def associate_iri_change_with_expression_change(summary_resting, summary_activated, title, fc_threshold, minRPKM, pc=0.001):
	"""
	Idea: we expect IR down would lead to expression up, IR up would lead to expression down 
	Approach: We find three groups of genes: iri_up, iri_down, iri_nochange. We then examine the expression changes of the three groups of genes. The expression changes are in terms of log2
	
	"""
	iri_up_genes = get_subset_genes_by_iri_changes(summary_resting, summary_activated, "up", fc_threshold, minRPKM, pc)
	iri_down_genes = get_subset_genes_by_iri_changes(summary_resting, summary_activated, "down", 1.0/fc_threshold, minRPKM, pc)
	iri_nochange_genes = get_subset_genes_by_iri_changes(summary_resting, summary_activated, "no-change", fc_threshold * 0.7, minRPKM, pc)
	
	expression_pc = 0.001
	expression_fc_for_iri_up_genes = get_expression_changes(summary_resting, summary_activated, iri_up_genes, expression_pc)#[(id, expression_fc)]
	expression_fc_for_iri_down_genes = get_expression_changes(summary_resting, summary_activated, iri_down_genes, expression_pc)#[(id, expression_fc)]
	expression_fc_for_iri_nochange_genes = get_expression_changes(summary_resting, summary_activated, iri_nochange_genes, expression_pc)#[(id, expression_fc)]
	
	mixed = []
	mixed.append([log(item[1],2) for item in expression_fc_for_iri_up_genes])
	mixed.append([log(item[1],2) for item in expression_fc_for_iri_nochange_genes])
	mixed.append([log(item[1],2) for item in expression_fc_for_iri_down_genes])
	
	#use R to plot
	title += "_associate_iri_change_with_expression_change"
	grdevices = importr('grDevices')
	grdevices.png(file = title + "-R.png", width=512, height=512)
	rlist = []
	for i in xrange(len(mixed)):
		rlist += r.list(robjects.FloatVector(mixed[i]))
		#rlist += robjects.FloatVector(mixed[i])
	#r.boxplot(rlist, col=r.c("gold", "darkgreen"), outline = r.FALSE)
	r.boxplot(rlist, names=r.c("iri_up", "iri_nochange", "iri_down"), col=r.c("dark orange", "dark green", "dark blue"))
	grdevices.dev_off()
	
	grdevices = importr('grDevices')
	grdevices.postscript(file = title + "-R.eps", width=512, height=512)
	r.boxplot(rlist, names=r.c("iri_up", "iri_nochange", "iri_down"), col=r.c("dark orange", "dark green", "dark blue"))
	grdevices.dev_off()
	
def find_candidate_IR_regulated_genes(summary_resting, summary_activated, fc_threshold, minRPKM, pc=0.001):
	"""
	Find genes with IR down and expression up, IR up with expression down
	Return: [ir_regulated_up_gene_ids], [ir_regulated_down_gene_ids]
	"""
	iri_up_gene_ids = get_subset_genes_by_iri_changes(summary_resting, summary_activated, "up", fc_threshold, minRPKM, pc)
	expresion_down_gene_ids = get_subset_genes_by_expression_changes(summary_resting, summary_activated, "down", 1.0/fc_threshold, minRPKM, pc)
	ir_regulated_up_gene_ids = list(set(iri_up_gene_ids) & set(expresion_down_gene_ids))
	
	iri_down_gene_ids = get_subset_genes_by_iri_changes(summary_resting, summary_activated, "down", 1.0/fc_threshold, minRPKM, pc)
	expresion_up_gene_ids = get_subset_genes_by_expression_changes(summary_resting, summary_activated, "up", fc_threshold, minRPKM, pc)
	ir_regulated_down_gene_ids = list(set(iri_down_gene_ids) & set(expresion_up_gene_ids))
	return ir_regulated_up_gene_ids, ir_regulated_down_gene_ids
		
def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--all_feature_pickle_file", action="store", type="string", dest="summary_file_names", help="names of summary pickle files of  RNASeq features", metavar="<file>")
	parser.add_option("-n", "--name", action="store", type="string", dest="title", help="name for the data set", metavar="<file>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print opt.species, " is not recognized, exiting";
		sys.exit(1);
	
	name_list, summary_files = get_names(opt.summary_file_names) #{name], {summary_name:file_name}
	my_summaries = {} #{summary_name:summary_object}
	for my_summary_name in summary_files.keys():
		my_summary_file_name = summary_files[my_summary_name]
		# load features for each summary
		if 	Utility_extended.fileExists(my_summary_file_name): #directly read the pickle file if it is available
			my_summary = pickle.load(open(my_summary_file_name, 'rb'))
			print "There are %d entrez IDs in %s" %(len(my_summary.keys()), my_summary_file_name)
			print "From the pickle there are %d features" % master.get_num_features(my_summary)
			print master.get_feature_names(my_summary)
			my_summaries[my_summary_name] = my_summary
		else:	
			print my_summary_file_name, " is not found";
			sys.exit(1)
		
	"""
	Calculate IRI
	"""
	iri_pc = 1 # read count pc
	for my_summary_name in my_summaries.keys():
		my_summary = my_summaries[my_summary_name]
		my_summary = copy.deepcopy(calculate_iri.calculate_expression_iri(my_summary, iri_pc))
		print "After adding iri. There are %d features in %s" % (master.get_num_features(my_summary), summary_files[my_summary_name])
		my_summaries[my_summary_name] = my_summary
		
	"""
	characterize IRI
	"""
	print "\nRemove genes with IRI>=1 and genes without intron"
	# Remove genes with likely misannotations (a misannotated id is an id with iri>1 in any of the libraries), remove genes with 0 introns
	max_iri = 1
	ids_misannotation_removed = removed_misannotated_ids_from_IRI(my_summaries, max_iri)
	print "\nAfter removing mis-annotations and genes with no introns, there are %d IDs left for downstream analysis" %len(ids_misannotation_removed)
	
	# Find expressed genes
	"""
	Three styles, 
	1) require that the genes are expressed in all states flag = my_summary.keys() 
	2) In a particular set of states, flag = subset of my_summary.keys() 
	3) In any one state,  flag =[-1] 
	"""
	expression_cutoff = 1
	flags = [name_list[0]]
	print "\nRemove genes with resting expression rpkm <1"
	expressed_ids = get_expressed_genes(my_summaries, expression_cutoff, flags)
	print "\nWith requirement of resting RPKM above ", expression_cutoff, ", ",  len(expressed_ids), " entrez IDs are left for downstream analysis"

	#Combine the requirements 
	print "\nRemove genes with resting expression rpkm <1 and with misannotation"
	subset_ids = list (set(ids_misannotation_removed) & set(expressed_ids))
	print "\nAfter removing mis-annotations,genes with no introns and genes not expressed, there are %d IDs left for downstream analysis" %(len(subset_ids))
	expressed_summaries = get_subset(my_summaries, subset_ids)
	
	
	print "\nCharacterize IRI between states for retained genes: boxplot"
	iri_pc = 0.001
	characterize_IRI(expressed_summaries, name_list, opt.title, iri_pc, None)	
	
	if len(name_list) == 2:
		fc_threshold = 2.0
		minRPKM = expression_cutoff
		pc = 0.001
		associate_iri_change_with_expression_change(expressed_summaries[name_list[0]], expressed_summaries[name_list[1]], opt.title, fc_threshold, minRPKM, pc)
	
		ir_regulated_up_gene_ids, ir_regulated_down_gene_ids = find_candidate_IR_regulated_genes(expressed_summaries[name_list[0]], expressed_summaries[name_list[1]], fc_threshold, minRPKM, pc=0.001)
	
if __name__ == "__main__":
	main(sys.argv)	
	
