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
import matplotlib

#import rpy2.robjects as robjects
#r = robjects.r
#from rpy2.robjects.packages import importr

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import Utility_extended
import GenomeData 
import SeparateByChrom
import UCSC_revised
import compare_two_libraries_on_genes
import gene_set_manipulation
import AnalyzeRNASeqPASeq

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")

# available from AnalyzeRNASeqPASeq
	# readfile(file)
	# write_entrez(entrez_genes, filename)
	# write_list(mylist, filename)
	# subfeature(entrez_genes, c)
	# subfeatures(entrez_genes, column_list)
	
	 
def subset_by_GeneSymbol(gene_symbols, entrez_genes, index=16):
	mysubset = {}
	for myid in entrez_genes.keys():
		my_gene_symbol = (entrez_genes[myid])[index]
		if my_gene_symbol in gene_symbols:
			mysubset[myid] = entrez_genes[myid]
	return mysubset

def getMainRefseqIDs(entrez_genes, outfile):
	mylist = [entrez_genes[myid][0] for myid in entrez_genes.keys()]
	out = open(outfile, 'w')
	for item in mylist:
		out.write(str(item) + '\n')
	out.close()
	return mylist

def getAllRefseqIDs(entrez_genes, outfile): 
	refseqs = []
	for myid in entrez_genes.keys():
		refseqs.extend(entrez_genes[myid][23].split(','))
	out = open(outfile, 'w')
	for item in refseqs:
		out.write(str(item) + '\n')
	out.close()
	return refseqs

def subset_by_values(entrez_genes, column_id, myrange):
	"""
	myrange: (start, end)
	Output an entrez_gene object
	"""
	mysubset = {}
	for myid in entrez_genes.keys():
		if inrange(float(entrez_genes[myid][column_id]), myrange) == 1:
			mysubset[myid] = entrez_genes[myid]
	return mysubset

def inrange(number, myrange):
	if number>=myrange[0] and number<myrange[1]:
		return 1
	else:
		return 0


def correlation(entrez_genes, c1, c2):
	mylist1 = []
	mylist2 = []
	for myid in entrez_genes.keys():
		mylist1.append(float((entrez_genes[myid])[c1]))
		mylist2.append(float((entrez_genes[myid])[c2]))
	print len(mylist1), len(mylist2)
	(pearson, spearman) = compare_two_libraries_on_genes.correlation(mylist1, mylist2, 0, 0)
	return (pearson, spearman)

def hist(entrez_genes, c1, title):
	"""
	need to make sure the column contains numbers
	"""
	mylist = []
	for myid in entrez_genes.keys():
		mylist.append(float((entrez_genes[myid])[c1]))
	myarray = numpy.array(mylist)
	hist = numpy.hist(myarray, bin=50, normed=True)
	plt.plot(hist, "bo", markersize=1.5);
	#plt.show()
	plt.savefig(title + ".png", format="png")
	plt.close()
	return hist		

def get_TSR(entrez_genes, downstream_extension=500, min_3UTR_length=50):
	"""
	TSR stands for the 3UTR shortening ratio. 
	"""
	resting = {}
	activated = {}
	for myid in entrez_genes.keys():
		threeUTR_length = float((entrez_genes[myid])[7])
		threeUTR_length -= downstream_extension
		resting_Length_Index = float((entrez_genes[myid])[1])
		activated_Length_Index = float((entrez_genes[myid])[2])
		
		#if threeUTR_length >= min_3UTR_length:
		resting_ratio = (threeUTR_length - resting_Length_Index)/threeUTR_length
		assert resting_ratio >= 0
		resting[myid] = resting_ratio
		
		activated_ratio = (threeUTR_length - activated_Length_Index)/threeUTR_length
		assert activated_ratio >= 0
		activated[myid] =  activated_ratio
	return (resting, activated)	


def get_intron_retention_index(entrez_genes, name="", stat = 0, pc = 0.000000000001):
	"""
	retains everything
	return dictionaries keyed by entrez_id, valued by iri
	if stat = 1, generate histogram
	"""
	resting_iri ={}
	activated_iri={}
	
	for myid in entrez_genes.keys():
		intron_length = float(entrez_genes[myid][22])
		if intron_length == 0:
			activated_iri[myid] = 0.0
			resting_iri[myid] = 0.0
		else:
			activated_iri[myid] = (float(entrez_genes[myid][21]) + pc)/(float(entrez_genes[myid][16]) + pc)
			resting_iri[myid] = (float(entrez_genes[myid][20]) + pc)/(float(entrez_genes[myid][15]) + pc)
			if float(entrez_genes[myid][20]) < 0.00000000001:
				resting_iri[myid] = 0
			if float(entrez_genes[myid][21]) < 0.00000000001:
				activated_iri[myid] = 0
	#font = {'family':'normal', 'weight':'bold', 'size':18}
	#matplotlib.rc('font', **font)
	if stat == 1 :
		log_resting_iri = numpy.array([log(i + pc, 2) for i in resting_iri.values()])
		log_activated_iri = numpy.array([log(i + pc, 2)  for i in activated_iri.values()])
		plt.hist(log_resting_iri, bins=100, color='r', normed=True, alpha=0.75,label="Resting" )
		plt.hist(log_activated_iri, bins=100, color='g', normed=True,  alpha=0.75, label="Activated")
		plt.title("IRI histogram")
		plt.xlabel("log(IRI)")
		plt.ylabel("Frequency")
		#plt.legend()
		plt.legend(loc = 'upper left')
		
		if name !="":
			plt.savefig( name + "activated_resting_iri_hist.eps", format="eps")
			plt.savefig( name + "activated_resting_iri_hist.png", format="png")
			
		
		resting_iri_limited = [i for i in resting_iri.values() if i>0 and i<1]
		print "Resting before log: mean", numpy.average(numpy.array(resting_iri_limited)), "    median:",  numpy.median(numpy.array(resting_iri_limited))
		activated_iri_limited = [i for i in activated_iri.values() if i>0 and i<1]
		print "Activated before log: mean", numpy.average(activated_iri_limited), "    median:", numpy.median(activated_iri_limited)	
		
		log_resting_iri = [log(i,2) for i in resting_iri_limited]
		log_activated_iri = [log(i,2) for i in activated_iri_limited]
		print "After log mean: resting ",  numpy.average(log_resting_iri), "    activated", numpy.average(log_activated_iri)
		print "After log median: resting ",  numpy.median(log_resting_iri), "    activated", numpy.median(log_activated_iri)
		print "After log std: resting ", numpy.std(log_resting_iri), "    activated", numpy.std(log_activated_iri)
		print scipy.stats.ttest_ind(log_resting_iri, log_activated_iri)
		print
		plt.clf()
		plt.hist(log_resting_iri, bins=50, color='r', normed=True, alpha=0.75,label="Resting" )
		plt.hist(log_activated_iri, bins=50, color='g', normed=True,  alpha=0.75, label="Activated")
		plt.title("IRI histogram")
		plt.xlabel("log(IRI)")
		plt.ylabel("Frequency")
		#plt.legend()
		plt.legend(loc = 'upper left')
		plt.savefig(name + "activated_resting_nonzero_lessthanone_iri_hist.eps", format="eps")
		plt.savefig(name + "activated_resting_nonzero_lessthanone_iri_hist.png", format="png")
		plt.close()
	return (resting_iri, activated_iri)

def get_subset_by_iri_ranking(entrez_genes, cell_type, ranking_type, number):
	"""
	rank the entrez_genes according to iri, and then pick a subset according to the ranking, like top 1000, bottom 1000, etc
	return an entrez_genes structure
	"""
	(resting_iri, activated_iri) = get_intron_retention_index(entrez_genes)
	if cell_type == "resting":
		mytuples = dic2list(resting_iri)
	elif cell_type == "activated":
		mytuples = dic2list(activated_iri)
	else:
		print "cell_type not recognized, must be resting or activated. "
		exit(1)
	sorted_tuples = sorted(mytuples, key=itemgetter(1))
	
	assert (len(mytuples) > 2 * number)
	if ranking_type == "top":
		myids = [item[0] for item in mytuples[-1*number:]]
	elif ranking_type == "bottom":
		myids = [item[0] for item in mytuples[0:number]]
	elif ranking_type == "middle":
		mid = len(mytuples)/2
		myrange = number/2
		assert (mid-myrange >= 0)
		assert (mid-myrange+number < len(mytuples))
		myids = [item[0] for item in mytuples[mid-myrange : mid-myrange+number]]
	else:
		print "ranking_type not recognized, must be top, middle or bottom. "
		exit(1)
	return subset(myids, entrez_genes)
	
def dic2list(dic):
	"""
	convert a dic of (key, value) to a list of tuples for sorting by value 
	"""
	mytuple = []
	for mykey in dic.keys():
		mytuple.append((mykey, dic[mykey]))
	return mytuple

def get_subset_genes_by_iri_changes(entrez_genes, change_type, fc_threshold, minRPKM, pc=0.001):
	"""
	change_type: "up", "down", "no-change"
	fc_threshold is activated_iri/resting_iri
	Only retain genes with shared-exon expression
	"""
	mysubset = {}
	for myid in entrez_genes.keys():
		
		resting_RPKM = float(entrez_genes[myid][15]) #shared exon
		activated_RPKM = float(entrez_genes[myid][16]) #shared exon
		if resting_RPKM > minRPKM and activated_RPKM > minRPKM:
			resting_iri = (float(entrez_genes[myid][20]) + pc)/float(resting_RPKM + pc)
			activated_iri = (float(entrez_genes[myid][21]) + pc)/float(activated_RPKM + pc)
			if change_type == "down":
				if activated_iri/resting_iri < fc_threshold:
					mysubset[myid] = entrez_genes[myid]
			elif change_type == "up":
				if activated_iri/resting_iri > fc_threshold:
					mysubset[myid] = entrez_genes[myid]
			elif change_type == "no-change":
				if activated_iri/resting_iri < fc_threshold and activated_iri/resting_iri > 1.0/fc_threshold:
					mysubset[myid] = entrez_genes[myid]
	return mysubset

	
def iri_RA_correlation(One, name):
	(resting_dic, activated_dic) = get_intron_retention_index(One)
	resting =  numpy.array(resting_dic.values()) 
	activated = numpy.array(activated_dic.values()) 
	print
	print "iri_RA_correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)
	
	# The plot clearly shows 4 clusters. 
	cluster0 = []
	assert (resting_dic.keys() == activated_dic.keys())
	for myid in resting_dic.keys():
		if resting_dic[myid]<10 and  activated_dic[myid]<10 :
			cluster0.append((myid,resting_dic[myid], activated_dic[myid]))
	print
	resting_cluster0 = [i[1] for i in cluster0]
	activate_cluster0 = [i[2] for i in cluster0]
	print "iri_RA_correlation cluster 0"
	print "Resting average: ", numpy.average(resting_cluster0), ".     Resting median: ", numpy.median(resting_cluster0)
	print "Activated average: ", numpy.average(activate_cluster0), ".     Activated median: ", numpy.median(activate_cluster0)
	scatterplot(resting_cluster0, activate_cluster0, name + "-cluster0")
	
def find_main_clusters_from_IRI(entrez_genes):
	"""
	scatter plot of resting vs activated produce 4 clusters. Three of the 4 have outrageous IRI in either of the cell type and the IRI is strongly dependent on pseudo count involved in defining IRI (intron RPKM + pc)/(exon RPKM + pc) . suggesting 1) mis-annotation and furthermore 2) for those genes, the exon read counts are 0, whereas the intron read counts are non-zero. For all meaningful analysis along the IRI line, we need to get rid of those genes. Perhaps it is useful to revisit the two differential IRI outlier clusters at a later stage. This module is for use of finding the IRI main cluster.
	
	The delination of the clusters are based on the IRI scatter plot:
	cluster 0 is the main cluster, independent of the pseudo count.
	cluster 1,2 IDs to be determined
	
	"""
	(resting_dic, activated_dic) = get_intron_retention_index(entrez_genes)
	assert (resting_dic.keys() == activated_dic.keys())
	cluster0 = []
	for myid in resting_dic.keys():
		if resting_dic[myid]<10 and  activated_dic[myid]<10 :
			cluster0.append(myid)
	return subset(cluster0, entrez_genes)
	
	

def get_resting_vs_activated_iri(entrez_genes):
	# Only retain genes with shared-exon expression and shared intron
	mylist = []
	(resting_iri, activated_iri) = get_intron_retention_index(entrez_genes)
	keys = (resting_iri.keys()).sort()
	for myid in keys:
		mylist.append((resting_iri[myid], activated_iri[myid]))
	return mylist

def get_irichanges(entrez_genes, minRPKM, pc=0.001):
	mysubset = {}
	for myid in entrez_genes.keys():
		resting_RPKM = float(entrez_genes[myid][15]) #shared exon
		activated_RPKM = float(entrez_genes[myid][16]) #shared exon
		if resting_RPKM > minRPKM or activated_RPKM > minRPKM:
			resting_iri = (float(entrez_genes[myid][20]) + pc)/float(resting_RPKM + pc)
			activated_iri = (float(entrez_genes[myid][21]) + pc)/float(activated_RPKM + pc)
			mysubset[myid] = activated_iri/resting_iri
	return mysubset

def irichange_vs_expressionchange(entrez_genes, minRPKM=0.5, pc = 0.0001):
	iri_changes = []
	expression_changes = []
	my_iri_changes = get_irichanges(entrez_genes, minRPKM, pc)
	my_expression_changes = get_expressionchanges(entrez_genes, minRPKM, pc)
	
	myidset = list( ( set(entrez_genes.keys()) & set(my_iri_changes.keys())) &  set(my_expression_changes.keys()))
	
	for myid in myidset:
		iri_changes.append(log(my_iri_changes[myid]))
		expression_changes.append(log(my_expression_changes[myid]))
		#iri_changes.append(my_iri_changes[myid])
		#expression_changes.append(my_expression_changes[myid])
	iri_changes_a = numpy.array(iri_changes)
	expression_changes_a = numpy.array(expression_changes)
	print "iri_changes(x) vs expression_changes(y) correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(iri_changes, expression_changes, 0, 0)
	#scatterplot(iri_changes_a, expression_changes_a, "iri_changes_vs_expression_changes")
	scatterplot(iri_changes_a, expression_changes_a, "iri_changes_vs_expression_changes", xscale='linear', yscale='linear')

def iri_vs_expressionchange(entrez_genes, minRPKM=0.5, pc = 0.0001):
	iri = []
	expression_changes = []
	(resting_iri, activated_iri) = get_intron_retention_index(entrez_genes, "", 0, pc)
	my_expression_changes = get_expressionchanges(entrez_genes, minRPKM, pc)
	
	myidset = list( ( set(entrez_genes.keys()) & set(resting_iri.keys())) &  set(my_expression_changes.keys()))
	
	for myid in myidset:
		iri.append(resting_iri[myid])
		expression_changes.append(log(my_expression_changes[myid]))
	iri_a = numpy.array(iri)
	expression_changes_a = numpy.array(expression_changes)
	print "iri(x) vs expression_changes(y) correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(iri, expression_changes, 0, 0)
	scatterplot(iri_a, expression_changes_a, "iri_vs_expression_changes")
	

def irichange_vs_3UTRIchange(entrez_genes, minRPKM=0.5, minreadcount=10, pc = 0.001):
	iri_changes = []
	threeUTR_changes = []
	my_iri_changes = get_irichanges(entrez_genes, minRPKM, pc)
	my_threeUTR_changes = get_threeUTRchanges(entrez_genes, minreadcount)
	
	myidset = list( ( set(entrez_genes.keys()) & set(my_iri_changes.keys())) &  set(my_threeUTR_changes.keys()))
	
	for myid in myidset:
		iri_changes.append(log(my_iri_changes[myid]))
		threeUTR_changes.append(my_threeUTR_changes[myid])
	iri_changes_a = numpy.array(iri_changes)
	threeUTR_changes_a = numpy.array(threeUTR_changes)
	print " log(iri_changes) vs 3UTR changes correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(iri_changes, threeUTR_changes, 0, 0)
	scatterplot(iri_changes_a, threeUTR_changes_a, "iri_changes_vs_threeUTR_changes", xscale='linear', yscale='linear')

def iri_vs_3UTRIchange(entrez_genes, minreadcount=10, pc = 0.001):
	iri = []
	threeUTR_changes = []
	(resting_iri, activated_iri) = get_intron_retention_index(entrez_genes, "", 0, pc)
	my_threeUTR_changes = get_threeUTRchanges(entrez_genes, minreadcount)
	
	myidset = list( ( set(entrez_genes.keys()) & set(resting_iri.keys())) &  set(my_threeUTR_changes.keys()))
	
	for myid in myidset:
		iri.append(resting_iri[myid])
		threeUTR_changes.append(my_threeUTR_changes[myid])
	iri_a = numpy.array(iri)
	threeUTR_changes_a = numpy.array(threeUTR_changes)
	print " log(resting iri) vs 3UTR changes correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(iri, threeUTR_changes, 0, 0)
	scatterplot(iri_a, threeUTR_changes_a, "resting_iri_vs_threeUTR_changes")

def iri_vs_irichange(entrez_genes, minRPKM=0.5, pc = 0.001):
	iri = []
	iri_changes = []
	(resting_iri, activated_iri) = get_intron_retention_index(entrez_genes, "", 0, pc)
	my_iri_changes = get_irichanges(entrez_genes, minRPKM, pc)
	
	myidset = list( ( set(entrez_genes.keys()) & set(resting_iri.keys())) &  set(my_iri_changes.keys()))
	
	for myid in myidset:
		iri.append(resting_iri[myid])
		iri_changes.append(my_iri_changes[myid])
			
	iri_a = numpy.array(iri)
	
	iri_changes_a = numpy.array(iri_changes)
	
	print "iri(x) vs iri_changes(y) correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(iri, iri_changes, 0, 0)
	#scatterplot(iri_changes_a, expression_changes_a, "iri_changes_vs_expression_changes")
	scatterplot(iri_a, iri_changes_a, "iri_vs_iri_changes")

def binning_by_expression(target_gene_exps, target_gene_names, all_genes_exp, number_of_bins):
	"""
	The module groups genes accord to expression, into number_of_bins bins
	output:a list of dictionaries, each of which is an entrezgenes object
	"""
	all_genes_feature = []
	for myid in all_genes_exp.keys():
		all_genes_feature.append((myid, float(all_genes_exp[myid][15]))) # shared exon RPKM
	all_genes_exps_grouped = Utility_extended.group2(all_genes_feature, 1, number_of_bins)
	
	target_genes_keys = target_gene_exps.keys()
	target_genes_binned = [] # a list of dictionaries, each of which is an entrez object
	
	for i in xrange(len(all_genes_exps_grouped)):
		this_group_keys = [item[0] for item in all_genes_exps_grouped[i]]
		this_group_target_keys = list( set(target_genes_keys) & set (this_group_keys) )
		target_genes_subset = subset(this_group_target_keys, all_genes_exp)
		target_genes_binned.append(target_genes_subset)
		#if len(this_group_target_keys) > 0:
			#my_target_filename = target_gene_names + "_group_" + str(i) + "_genes"
			#write_entrez(target_genes_subset, my_target_filename + ".dat")
		#this_group_remaining_keys = list( set (this_group_keys) - set(this_group_target_keys) )
		#control_genes_subset = subset(this_group_remaining_keys, all_genes_exp)
		#control_genes_binned.append(control_genes_subset)
		#if len(this_group_remaining_keys) > 0:
			#my_control_genes_filename = remaining_gene_names + "_group_" + str(i) + "_genes"
			#write_entrez(control_genes_subset, my_control_genes_filename+".dat")
	return target_genes_binned



def compare_group_expressions(target_genes_binned, control_genes_binned, pc = 0.09):
	
	assert (len(target_genes_binned) == len(control_genes_binned))
	target_gene_resting = []
	control_genes_resting = []
	target_gene_activated = []
	control_genes_activated = []
	for i in xrange(len(target_genes_binned)):
		target_gene_resting.append([log(float((target_genes_binned[i])[myid][15]) + pc,2) for myid in target_genes_binned[i].keys()])
		control_genes_resting.append([log(float((control_genes_binned[i])[myid][15]) + pc,2) for myid in control_genes_binned[i].keys()])
		target_gene_activated.append([log(float((target_genes_binned[i])[myid][16]) + pc,2) for myid in target_genes_binned[i].keys()])
		control_genes_activated.append([log(float((control_genes_binned[i])[myid][16]) + pc,2) for myid in control_genes_binned[i].keys()])
	plt.figure()
	comparative_boxplot(target_gene_resting, control_genes_resting, "resting_expression_target_vs_remaining")
	comparative_boxplot(target_gene_activated, control_genes_activated, "activated_expression_target_vs_remaining")
	comparative_boxplot(target_gene_resting, target_gene_activated, "target_expression_resting_vs_activated")
	comparative_boxplot(control_genes_resting, control_genes_activated, "remaining_expression_resting_vs_activated")

def comparative_boxplot(target_genes_binned, control_genes_binned, title):
	"""
	both are list of lists of numbers
	"""
	print title
	assert (len(target_genes_binned) == len(control_genes_binned))
	
	#plt.clf()
	mixed = []
	for i in xrange(len(target_genes_binned)):
		mixed.append(target_genes_binned[i])
		mixed.append(control_genes_binned[i])
		if len(target_genes_binned[i]) > 0:
			print i,"th group: "
			print "		Number of genes: ", len(target_genes_binned[i]), len(control_genes_binned[i])
			print "		p-value associated with difference in distribution using k-s test: ", scipy.stats.ks_2samp(target_genes_binned[i], control_genes_binned[i])[1]	
			print "		p-value associated with ranking of the two sets of observations using MWU test: ", scipy.stats.mannwhitneyu(control_genes_binned[i], target_genes_binned[i])[1]	
	#myboxplot = plt.boxplot(mixed)
	#plt.setp(myboxplot['medians'], color = 'green')
	#plt.setp(myboxplot['boxes'], color = 'black')
	#plt.setp(myboxplot['whiskers'], color = 'black')
	#plt.setp(myboxplot['fliers'], color = 'green')
	#plt.savefig(title + ".png", format="png")
	#use R to plot
	grdevices = importr('grDevices')
	grdevices.png(file = title + "-R.png", width=512, height=512)
	rlist = []
	for i in xrange(len(mixed)):
		rlist += r.list(robjects.FloatVector(mixed[i]))
		#rlist += robjects.FloatVector(mixed[i])
	r.boxplot(rlist, col=r.c("gold", "darkgreen"))
	grdevices.dev_off()
	
	return mixed

def compare_feature_on_binned_genes(target_genes_binned, control_genes_binned, feature, title):
	"""
	feature is a dic keyed by entrez id
	"""
	assert (len(target_genes_binned) == len(control_genes_binned))
	feature_on_target = []
	feature_on_control = []
	for i in xrange(len(target_genes_binned)):
		target_genes = target_genes_binned[i]
		mykeys = list ( set(target_genes.keys()) & set(feature.keys()))
		feature_on_target.append([ feature[mykey] for mykey in mykeys])
		
		control_genes = control_genes_binned[i]
		mykeys = list ( set(control_genes.keys()) & set(feature.keys()))
		feature_on_control.append([ feature[mykey] for mykey in mykeys])

	mixed = comparative_boxplot(feature_on_target, feature_on_control,  title)	
	return mixed

def compare_feature_on_sets_of_binned_genes(sets_of_binned_genes, feature, title):
	"""
	sets_of_binned_genes: a tuple of list, each list represent a specific subset of genes binned by expression.  
	feature is a dic keyed by entrez id
	"""
	feature_on_sets_of_binned_genes = []
	for i in xrange(len(sets_of_binned_genes)):
		assert (len(sets_of_binned_genes[0]) == len(sets_of_binned_genes[i])) #len should be number of bins
	for i in xrange(len(sets_of_binned_genes)):
		feature_on_sets_of_binned_genes.append([])
		for j in xrange(len(sets_of_binned_genes[i])):
			target_genes = (sets_of_binned_genes[i])[j]
			mykeys = list ( set(target_genes.keys()) & set(feature.keys()))
			feature_on_sets_of_binned_genes[i].append([ feature[mykey] for mykey in mykeys])
	mixed = comparative_multi_boxplot(feature_on_sets_of_binned_genes,  title)	
	return mixed
	
def comparative_multi_boxplot (feature_on_sets_of_binned_genes, title):
	"""
	feature_on_sets_of_binned_genes: a tuple, with element being list, whose element is a list of numbers
	"""
	mixed = []
	length = len(feature_on_sets_of_binned_genes)
	number_of_bins = len(feature_sets_of_binned_genes[0])
	for i in xrange(len(feature_sets_of_binned_genes)):
		assert ((number_of_bins) == len(feature_sets_of_binned_genes[i])) #len should be number of bins
	for j in xrange(number_of_bins):
		for i in xrange(len(feature_sets_of_binned_genes)):
			mixed.append((feature_sets_of_binned_genes[i])[j])
	grdevices = importr('grDevices')
	grdevices.png(file = title + "-R.png", width=512, height=512)
	rlist = []
	mycolor = r.c("gold", "darkgreen", "light blue", "pink", "brown", "dark red")
	for i in xrange(len(mixed)):
		rlist += r.list(robjects.FloatVector(mixed[i]))
	r.boxplot(rlist, col=mycolor[0:length+1])
	grdevices.dev_off()
	return mixed

def compare_histonemark_on_sets_of_binned_genes(sets_of_binned_genes, histonemod, title, histonemod_ind=11, pc=0.009):
	"""
	The module tries to dissect the difference between sets of genes and control genes from histonemod.	
	sets_of_binned_genes are lists of entrez_genes objects
	"""
	histonmod_easy={}
	for mykey in histonemod.keys():
		histonmod_easy[mykey] = log(float((histonemod[mykey])[histonemod_ind]) + pc ,2)
	mixed = compare_feature_on_sets_of_binned_genes(sets_of_binned_genes, histonmod_easy, title)
	return mixed 

def compare_expression_on_sets_of_binned_genes(sets_of_binned_genes, pc = 0.09):
	resting_expressions = []
	activated_expressions = []
	number_of_bins = len(sets_of_binned_genes[0])
	for i in xrange(len(sets_of_binned_genes)):
		assert ((number_of_bins) == len(sets_of_binned_genes[i])) #len should be number of bins
	for i in xrange(len(sets_of_binned_genes)):
		resting_expressions.append([])
		activated_expressions.append([])
		for j in xrange(number_of_bins):
			resting_expressions[i].append([log(float((sets_of_binned_genes[i][j])[myid][15]) + pc,2) for myid in (sets_of_binned_genes[i])[j].keys()]) 
			activated_expressions[i].append([log(float((sets_of_binned_genes[i][j])[myid][16]) + pc,2) for myid in (sets_of_binned_genes[i])[j].keys()])
	comparative_multi_boxplot(resting_expressions, "resting_expression")
	comparative_multi_boxplot(activated_expressions, "activated_expression")
	for i in xrange(len(sets_of_binned_genes)):
		comparative_boxplot(resting_expressions[i], activated_expressions[i], "resting_vs_activated_" + str(i))

def get_mark_values(genes, histonemod, histonemod_ind=11):
	"""
	histonemod_ind 11 transcript
	Return a list of tuples (ID, mark_value)
	"""
	mark_values = []
	keys = list ( set(genes.keys()) & set(histonemod.keys()))
	for ID in keys:
		mark_values.append((ID, float(histonemod[ID][histonemod_ind])))
	return mark_values

def compare_histonemark_on_binned_genes(target_genes_binned, control_genes_binned, histonemod, title, histonemod_ind=11, pc=0.009):
	"""
	The module tries to dissect the difference between target genes and control genes from histonemod.	
	target_genes_binned, control_genes_binned are lists of entrez_genes objects
	"""
	histonmod_easy={}
	for mykey in histonemod.keys():
		histonmod_easy[mykey] = log(float((histonemod[mykey])[histonemod_ind]) + pc ,2)
	mixed = compare_feature_on_binned_genes(target_genes_binned, control_genes_binned, histonmod_easy, title)
	return mixed 


def get_epmark_iri(histmod, pc=0.00001):
	iri = {} # iri keyed by entrez id and valued by iri
	entrez_ids = histmod.keys()
	excluded = 0
	if len(entrez_ids) > 0:
		#test_id = entrez_ids[:10]
		for myid in entrez_ids:
			intron_length = float(histmod[myid][7])
			if intron_length == 0:
				iri[myid] = 0.0
			elif float(histmod[myid][8]) < 0.0000000001:
				iri[myid] = 0
			else:
				if float(histmod[myid][5]) == 0: # exon read count is zero
					excluded += 1
				else:
					iri[myid] = (float(histmod[myid][8]) + pc)/(float(histmod[myid][5]) + pc)
			#if myid in test_id:
				#print  histmod[myid]
				#print intron_length
				#print histmod[myid][8], histmod[myid][5], iri[myid]
				#print
		iri_list = numpy.array(iri.values())
		
		print "the number of genes with zero exons read count: ", excluded 
		print "IRI Average is: ", numpy.average(iri_list)
		print "IRI Median is: ", numpy.median(iri_list)
	return iri

def compare_marks_iri_on_binned_genes(target_genes_binned, control_genes_binned, histmod, histmod_name, pc= 0.00001):
	
	iri = get_epmark_iri(histmod, pc)
	for myid in iri.keys():
		iri[myid] = log(iri[myid]+pc, 2)
	title = histmod_name + "_IRI_target_vs_control_boxplot"
	mixed = compare_feature_on_binned_genes(target_genes_binned, control_genes_binned, iri, title)
	return mixed

def compare_iri_histone_vs_expression(entrez_genes, histmod, title):
	pc = 0.00001
	resting_exp_iri = (get_intron_retention_index(entrez_genes, "", 0, pc))[0] # dic keyed by entrez id, valued by iri
	resting_hist_iri = get_epmark_iri(histmod, pc)
	myid_set = list(set(resting_exp_iri.keys()) & set(resting_hist_iri.keys()))
	resting_exp_iri_list = []
	resting_hist_iri_list = []
	for myid in myid_set:
		resting_exp_iri_list.append(resting_exp_iri[myid])
		resting_hist_iri_list.append(resting_hist_iri[myid])
	compare_two_libraries_on_genes.correlation(resting_exp_iri_list, resting_hist_iri_list, 0, 0)
	scatterplot (resting_exp_iri_list, resting_hist_iri_list, title)	

def compare_PAReadCount_on_binned_genes(target_genes_binned, control_genes_binned, entrez_genes, celltype, pc = 0.1):
	feature = {} 
	if celltype == "resting":
		for myid in entrez_genes.keys():
			feature[myid] = log(float((entrez_genes[myid])[5]) + pc,2) # resting 3UTR PA read count 
	elif celltype == "activated":
		for myid in entrez_genes.keys():
			feature[myid] = log(float((entrez_genes[myid])[6])+ pc, 2) # activated 3UTR PA read count 
	else:
		print "cell-type has to be resting or activated"
		exit(1)
	title = celltype + "_PA_3UTR_read_count" + "_target_vs_control_boxplot"
	mixed = compare_feature_on_binned_genes(target_genes_binned, control_genes_binned, feature, title)

def compare_3UTRshortenningenrichment_on_binned_genes(target_genes_binned, control_genes_binned, cutoff, title):
	"""
	cutoff is the minimum distance for defined shortening 
	"""
	shortened_target_fractions = []
	shortened_control_fractions = []
	for i in xrange(len(target_genes_binned)):
		target_genes = target_genes_binned[i]
		control_genes = control_genes_binned[i]
		#print "\t number of target genes: ", len(target_genes)
		#print "\t number of control genes: ", len(control_genes)	

		total = float(len(target_genes.values()))
		TUTRshortened_target_genes = get_3UTRshortened_genes(target_genes, 0, cutoff)
		shortened_target_fraction.append([len(TUTRshortened_target_genes.values())/total])

		total = float(len(control_genes.values()))
		TUTRshortened_control_genes = get_3UTRshortened_genes(control_genes, 0, cutoff)
		shortened_control_fraction.append([len(TUTRshortened_control_genes.values())/total])
	mixed = comparative_boxplot(shortened_target_fractions, shortened_control_fractions, title)	
	return mixed

def get_expressionchanges(entrez_genes, minRPKM, pc = 0.001):
	myentrez={}
	for myid in entrez_genes.keys():
		resting_RPKM = float(entrez_genes[myid][15])
		activated_RPKM = float(entrez_genes[myid][16])
		fc = (activated_RPKM + pc)/(resting_RPKM + pc)
		if activated_RPKM >= minRPKM or resting_RPKM >= minRPKM:
			myentrez[myid] = fc
	return myentrez

def get_expressed_genes(entrez_genes, state,  minRPKM):
	myentrez={}
	for myid in entrez_genes.keys():
		if state == "resting":
			resting_RPKM = float(entrez_genes[myid][15])
			if resting_RPKM >= minRPKM:
				myentrez[myid] = entrez_genes[myid]
		elif state == "activated":
			activated_RPKM = float(entrez_genes[myid][16])
			if activated_RPKM >= minRPKM:
				myentrez[myid] = entrez_genes[myid]
		else:
			print "State can only be resting or activated"
			exit(1)
	return myentrez

def get_expression_changed_genes(entrez_genes, fc_threshold, minRPKM, pc = 0.001):
	"""
	activated vs resting
	"""
	myentrez={}
	for myid in entrez_genes.keys():
		resting_RPKM = float(entrez_genes[myid][15])
		activated_RPKM = float(entrez_genes[myid][16])
		fc = (activated_RPKM + pc)/(resting_RPKM + pc)
		if fc_threshold >= 1:
			if activated_RPKM >= minRPKM and fc >= fc_threshold:
				myentrez[myid] = entrez_genes[myid]
		elif fc_threshold < 1:
			if resting_RPKM >= minRPKM and fc <= fc_threshold:
				myentrez[myid] = entrez_genes[myid]
	return myentrez

def get_threeUTRchanges(entrez_genes, minreadcount):
	"""
	From PASeq results
	"""
	mysubset = {}
	
	for myid in entrez_genes.keys():
		resting_length_index = float((entrez_genes[myid])[1])
		resting_PA_rc = float((entrez_genes[myid])[5])
		activated_length_index = float((entrez_genes[myid])[2])
		activated_PA_rc = float((entrez_genes[myid])[6])
		if min(resting_PA_rc, activated_PA_rc)  > minreadcount:
			threeUTR_length = float((entrez_genes[myid])[7])
			#mysubset[myid] = (activated_length_index - resting_length_index)/float(threeUTR_length)
			mysubset[myid] = (activated_length_index - resting_length_index)
	return mysubset
	
def get_3UTRshortened_genes(entrez_genes, minreadcount, min_distance):
	"""
	For PASeq results
	"""
	mysubset = {}
	for myid in entrez_genes.keys():
		resting_length_index = float((entrez_genes[myid])[1])
		resting_PA_rc = float((entrez_genes[myid])[5])
		activated_length_index = float((entrez_genes[myid])[2])
		activated_PA_rc = float((entrez_genes[myid])[6])
		if min(resting_PA_rc, activated_PA_rc)  > minreadcount and  (activated_length_index - resting_length_index) > min_distance:
			mysubset[myid] = entrez_genes[myid]
	return mysubset

def get_3UTRshortened_PMInochange_genes(entrez_genes, minreadcount, min_distance):
	tUTRshortened = get_3UTRshortened_genes(entrez_genes, minreadcount, min_distance)
	mysubset = {}
	for myid in entrez_genes.keys():
		resting_PMI= float((entrez_genes[myid])[3])
		activated_PMI = float((entrez_genes[myid])[4])
		if (activated_PMI-resting_PMI)/(activated_PMI+0.1) < 1:
			mysubset[myid] = entrez_genes[myid]
	myset = list(set(tUTRshortened.keys()) & set(mysubset.keys()))
	print len(myset)
	my_entrez_subset = subset(myset, entrez_genes)
	write_entrez(my_entrez_subset, "3UTRshortened_PMInochange_genes.txt")
	return my_entrez_subset

def expressionchange_vs_3UTRIchange (entrez_genes, minRPKM=0.5, pc = 0.001, minreadcount=50):
	tUTRI_changes = []
	expression_changes = []
	my_tUTRI_changes = get_threeUTRchanges(entrez_genes, minreadcount)
	my_expression_changes = get_expressionchanges(entrez_genes, minRPKM, pc)
	myidset = list( ( set(entrez_genes.keys()) & set(my_tUTRI_changes.keys())) &  set(my_expression_changes.keys()))
	for myid in myidset:
		tUTRI_changes.append(my_tUTRI_changes[myid])
		expression_changes.append(log(my_expression_changes[myid], 2))
		#tUTRI_changes.append(my_tUTRI_changes[myid])
		#expression_changes.append(my_expression_changes[myid])
	tUTRI_changes_a = numpy.array(tUTRI_changes)
	expression_changes_a = numpy.array(expression_changes)
	print "tUTRI_changes(x) vs expression_changes(y) correlations (pearson, spearman): ",  compare_two_libraries_on_genes.correlation(tUTRI_changes, expression_changes, 0, 0)
	scatterplot(expression_changes_a, tUTRI_changes_a, "expression_changes_vs_3UTR_fraction_changes", xscale='linear', yscale='linear')

def scatterplot(a, b, title, xscale='log', yscale='log'):
	plt.plot(a, b, "bo", markersize = 1.5, alpha = 1);
	ax = plt.gca();
	ax.set_xscale(xscale)
	ax.set_yscale(yscale)
	#ax.set_aspect(1.) 
	ax.grid (color='gray', linestyle='dashed')
	#plt.savefig(title + ".eps", format="eps")
	plt.savefig(title + ".png", format="png")
	plt.close()
	
def merged_exons_RPKM_RA_correlation(One, name):
	resting =  numpy.array(subfeature(One, 10)) #merged exons
	activated = numpy.array(subfeature(One, 11)) #merged exons
	print "Merge exons expression RPKM RA correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)

def shared_exons_RPKM_RA_correlation(One, name):
	resting =  numpy.array(subfeature(One, 15)) #merged exons
	activated = numpy.array(subfeature(One, 16)) #merged exons
	print
	print "Shared exons expression RPKM RA correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)

def ThreeUTRI_RA_correlation(One, name):
	"""
	Resting vs Activated
	"""
	resting = numpy.array(subfeature(One, 1))
	activated = numpy.array(subfeature(One, 2))
	print
	print "ThreeUTRI_RA_correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)

def PMI_RA_correlation(One, name):
	"""
	RA: Resting vs activated
	"""
	resting = numpy.array(subfeature(One, 3))
	activated = numpy.array(subfeature(One, 4))
	print
	print "PMI_RA_correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)

def TSR_RA_correlation(One, name, downstream_extension=500, min_3UTR_length=50):
	(resting_dic, activated_dic) = get_TSR(One, downstream_extension, min_3UTR_length)
	resting = numpy.array(resting_dic.values())
	activated = numpy.array(activated_dic.values())
	print
	print "TSR_RA_correlation"
	print "Resting average: ", numpy.average(resting), ".     Resting median: ", numpy.median(resting)
	print "Activated average: ", numpy.average(activated), ".     Activated median: ", numpy.median(activated)
	scatterplot(resting, activated, name)



def ThreeUTRI_PMI_correlation(One, name):
	"""
	What could result in a linear correlation?
	~2 sites per gene?
	"""
	#print 
	#print "ThreeUTRI_PMI_correlation"
	resting_3UTRI = numpy.array(subfeature(One, 1))
	activated_3UTRI = numpy.array(subfeature(One, 2))
	resting_PMI = numpy.array(subfeature(One, 3))
	activated_PMI = numpy.array(subfeature(One, 4))
	scatterplot(resting_3UTRI, resting_PMI, "Resting-" + name)
	scatterplot(activated_3UTRI, activated_PMI, "Activated-" + name)

def ThreeUTRI_IRI_correlation(One, name):
	resting_3UTRI = numpy.array(subfeature(One, 1))
	activated_3UTRI = numpy.array(subfeature(One, 2))
	(resting_dic, activated_dic) = get_intron_retention_index(One)
	resting_IRI =  numpy.array(resting_dic.values()) 
	activated_IRI = numpy.array(activated_dic.values()) 
	scatterplot(resting_3UTRI, resting_IRI, "Resting-" + name)
	scatterplot(activated_3UTRI, activated_IRI, "Activated-" + name)

def TSR_IRI_correlation(One, name):
	(resting_dic, activated_dic) = get_TSR(One)
	resting_TSR = numpy.array(resting_dic.values())
	activated_TSR = numpy.array(activated_dic.values())
	(resting_dic, activated_dic) = get_intron_retention_index(One)
	resting_IRI =  numpy.array(resting_dic.values()) 
	activated_IRI = numpy.array(activated_dic.values()) 
	scatterplot(resting_TSR, resting_IRI, "Resting-" + name)
	scatterplot(activated_TSR, activated_IRI, "Activated-" + name)

def scatter_correlation(entrez_genes, geneset_name):
	"""
	Correlation in on cell type
	cell_type = "resting" or "activated"
	"""
	# Calculate correlation of PA vs RNA (merged exons) in resting 
	print "Correlation of PA vs RNA (merged exons) in resting: ", correlation(entrez_genes, 10, 5)
	print "Correlation of PA vs RNA (merged exons) in activated: ", correlation(entrez_genes, 11, 6)
	
	# Calculate correlation of PA vs RNA (shared exons) in resting 
	print "Correlation of PA vs RNA (shared exons) in resting: ", correlation(entrez_genes, 15, 5)
	print "Correlation of PA vs RNA (shared exons) in activated: ", correlation(entrez_genes, 16, 6)
	
	resting_RPKM =  numpy.array(subfeature(entrez_genes, 10)) #merged exons
	activated_RPKM = numpy.array(subfeature(entrez_genes, 11)) #merged exons
	resting_3UTRI = numpy.array(subfeature(entrez_genes, 1))
	activated_3UTRI = numpy.array(subfeature(entrez_genes, 2))
	scatterplot(resting_RPKM, resting_3UTRI, "Resting_RNA_3UTRI_" + geneset_name)
	scatterplot(activated_RPKM, activated_3UTRI, "Activated_RNA_3UTRI_" + geneset_name)

def scatter_RA_correlations(One, geneset_name):
	merged_exons_RPKM_RA_correlation(One, "Merged_exons_RPKM_" + geneset_name)
	#shared_exons_RPKM_RA_correlation(One, "Shared_exons_RPKM_" + geneset_name)
	iri_RA_correlation(One, "Iri_" + geneset_name)
	ThreeUTRI_RA_correlation(One, "ThreeUTRI_" + geneset_name)
	PMI_RA_correlation(One, "PMI_" + geneset_name)
	TSR_RA_correlation(One, "TSI_" + geneset_name, 500, 50)
	ThreeUTRI_PMI_correlation(One, "3UTR_vs_PMI_" + geneset_name)
	ThreeUTRI_IRI_correlation(One, "3UTR_vs_IRI_" + geneset_name)
	TSR_IRI_correlation(One,"TSR_vs_IRI_" + geneset_name)

def global_survey(entrez_genes, geneset_name):
	scatter_correlation(entrez_genes, geneset_name)
	scatter_RA_correlations(entrez_genes, geneset_name)

def IdentifyTargetControl_1(entrez_genes, number_of_bins):
	"""
	target: IRI high and down genes
	target_genes_name = "IRI_high_and_down_genes"	
	control: rest
	"""
	iri_down_genes = get_subset_genes_by_iri_changes(entrez_genes, "down", 0.5, 0.15, 0.0000001)
	print "Number of iri_down_genes: ", len(iri_down_genes.keys())
	write_entrez(iri_down_genes, "iri_down_genes.dat")
	#iri_up_genes = get_subset_genes_by_iri_changes(entrez_genes, 'up', 2, 0.15, 0.0000001)
	#print "iri_up_genes: ",  len(iri_up_genes.keys())
	# compare the epimarks from target genes and the control genee. 
	
	resting_iri = get_intron_retention_index(entrez_genes)[0]
	target_genes = {}
	for myid in resting_iri.keys():
		# These choices results in 1600 genes, whose expression has a similar spread as those of all genes.
		if resting_iri[myid] < 1 and resting_iri[myid] >= exp(-2):
			target_genes[myid] = expressed_one[myid] #IRI high genes
	print "The number of IRI high genes: ", len(target_genes.keys())
	write_entrez(target_genes, "IRI_high_genes.txt")		
	
	target_genes_name = "IRI_high_and_down_genes"
	#The target genes are IRI high and down genes, and control genes are the rest.  
	myid_set = list(set(target_genes.keys()) & set(iri_down_genes.keys()))
	print "Number of IRI_high_and_down_genes: ", len(myid_set)
	target_genes_set = subset(myid_set, entrez_genes)
	write_entrez(target_genes_set, "IRI_high_and_down_genes.dat")
	#irichange_vs_expressionchange(target_genes_set, 0.15)
	
	target_genes_binned = binning_by_expression(target_genes_set, target_genes_name, entrez_genes, number_of_bins)
	
	control_genes_ids = list( set(entrez_genes.keys()) - set(target_genes_set.keys()))
	control_genes_name = "Not_" + target_genes_name
	control_genes_binned = binning_by_expression(control_genes_set, control_genes_name, entrez_genes, number_of_bins)	
	
	return (target_genes_binned, control_genes_binned)



def IdentifyTargetControl_2(entrez_genes, number_of_bins):
	"""
	target: IRI high and down genes
	target_genes_name = "IRI_high_and_down_genes"	
	control: rest
	"""
	iri_top1000 = get_subset_by_iri_ranking(entrez_genes, "resting", "top", 1000)
	iri_bottom1000 = get_subset_by_iri_ranking(entrez_genes, "resting", "bottom", 1000)
	
	

#def top_bottom_comparison(entrez_genes, histmod, number_of_bins, histonemod_ind=11, pc=0.001):

	##compare_group_expressions(target_genes_binned, control_genes_binned)
	##title = "resting_" + libName1 + "_target_vs_control"
	##compare_histonemark_on_binned_genes(target_genes_binned, control_genes_binned, HistoneMod, title, column_id)
	##compare_marks_iri_on_binned_genes(target_genes_binned, control_genes_binned, HistoneMod, libName1, 0.01)
	##compare_PAReadCount_on_binned_genes(target_genes_binned, control_genes_binned, expressed_one, "resting")


def Analyzegenesets(sets_of_binned_genes):
	"""
	sets_of_binned_genes: tuple of gene lists
	"""
	



def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--combined_result_file", action="store", type="string", dest="one", help="combined result file", metavar="<file>")
	parser.add_option("-j", "--histone_modification_file", action="store", type="string", dest="HistoneMod", help="HistoneMod at resting", metavar="<file>")
	parser.add_option("-c", "--histone_modification_location", action="store", type="string", dest="histone_modification_location", help="transcript, shared_exon, merged_exon, intron", metavar="<str>")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
		
	if Utility_extended.fileExists(opt.HistoneMod):
		HistoneMod = readfile(opt.HistoneMod)
		# print len(HistoneMod.keys())
		# Read the peaks, which is assumed to have the pseudo ucsc format
		libName1 = (opt.HistoneMod).split('/')[-1]
		libName1 = libName1.split('.')[0]
		libName2 = libName1
		
		if opt.histone_modification_location == "transcript":
			column_id = 11
		elif opt.histone_modification_location == "shared_exon":
			column_id = 5
		elif opt.histone_modification_location == "merged_exon":
			column_id = 2
		elif opt.histone_modification_location == "intron":
			column_id = 8
		else:
			print "wrong choice for ", opt.histone_modification_location
			exit(1)
		
	else:
		print opt.HistoneMod, " is not found";
		sys.exit(1)
	
	
	TA = {}
	TA[0] = "Entrez ID" 	 
	TA[0] = "Main Refseq ID"	
	TA[1] = "Resting Length Index"	
	TA[2] =	"Activated Length Index"	
	TA[3] =	"Resting PA Multiplicity Index"	
	TA[4] =	"Activated PA Multiplicity Index"	
	TA[5] = "Resting Peaks Read Count"	
	TA[6] = "Activated Peaks Read Count"	
	TA[7] =	"3UTR length"	 
	TA[8] = "Resting Merged Exon Read Count"	 
	TA[9] = "Activated Merged Exon Read Count"	
	TA[10] = "Resting Merged Exon RPKM"	
	TA[11] = "Activated Merged Exon RPKM"	
	TA[12] = "Merged Exon Length"	
	TA[13] = "Resting Shared Exon Read Count"	 
	TA[14] = "Activated Shared Exon Read Count" 	
	TA[15] = "Resting Shared Exon RPKM"	
	TA[16] = "Activated Shared Exon RPKM"	
	TA[17] = "Shared Exon Length"	
	TA[18] = "Resting Shared Intron Read Count" 	
	TA[19] = "Activated Shared Intron Read Count" 	
	TA[20] = "Resting Shared Intron RPKM"	
	TA[21] = "Activated Shared Intron RPKM"
	TA[22] = "Share Intron Length"
	TA[23] = "RefSeq IDs"	 
	TA[24] = "Gene Symbols"
	
	
	# choose the starting geneset
	One = readfile(opt.one)
	lib_RNASeq_name = (opt.one).split('/')[-1]
	lib_RNASeq_name = lib_RNASeq_name.split('.')[0]
	print "There are ", len(One.keys()), " entrez IDs from the outset"
	Misannotation_removed_One = find_main_clusters_from_IRI(One)
	write_entrez(Misannotation_removed_One, opt.one + "-misannotation-removed")
	print "After removing mis-annotations, there are ", len(Misannotation_removed_One.keys()), " entrez IDs are left for downstream analysis"
	cutoff = 0.15
	expressed_one = get_expressed_genes(Misannotation_removed_One, 'resting', cutoff)
	print "After expression requirement of above ", cutoff, ", ",  len(expressed_one.keys()), " entrez IDs are left for downstream analysis"

	# 
	libName1 = libName1 +"_on_" + opt.histone_modification_location
	title = "resting_expiri_" + libName1 + "_iri_all_genes"
	#compare_iri_histone_vs_expression(Misannotation_removed_One, HistoneMod, title)
	
	(resting_iri, activated_iri) = get_intron_retention_index(expressed_one, lib_RNASeq_name +"_expressed_",1)
	write_entrez(resting_iri, "Resting_IRI_expressed_genes.txt")
	write_entrez(activated_iri, "Activated_IRI_expressed_genes.txt")
	exit(1)
	
	#title = "resting_expiri_histiri_target_genes"
	#compare_iri_histone_vs_expression(target_genes_set, HistoneMod,title)
	

	
	#compare_group_expressions(target_genes_binned, control_genes_binned)
	#title = "resting_" + libName1 + "_target_vs_control"
	#compare_histonemark_on_binned_genes(target_genes_binned, control_genes_binned, HistoneMod, title, column_id)
	#compare_marks_iri_on_binned_genes(target_genes_binned, control_genes_binned, HistoneMod, libName2, 0.01)
	#compare_PAReadCount_on_binned_genes(target_genes_binned, control_genes_binned, expressed_one, "resting")

	
	
if __name__ == "__main__":
	main(sys.argv)