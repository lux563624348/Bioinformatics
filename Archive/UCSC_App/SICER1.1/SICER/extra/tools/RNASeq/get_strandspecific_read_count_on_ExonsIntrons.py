#!/usr/bin/env python

"""

Calculate strand-specific read counts on intronic, exonic regions with strand specific reads (eg, strand-specific RNA Seq)
Output:
1) .dat file that contains 
outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"

2) pickle files

2.1) 	summary: {entrezID:{attribute:value}}
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

2.2	all_reads_on_shared_exons = {} # {entrezID:[((start, end), read_count)]}
	rc on individual exons

2.3	all_reads_on_shared_introns = {} # {entrezID:[((start, end), read_count)]}
	rc on individual introns
	
2.4	all_reads_on_merged_transcripts = {} #{entrezID:[((start, end), read_count)]}
	rc on the merged transcript, expected to have only one element for each entrezID		

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import copy
import numpy
import matplotlib.pyplot as plt
import matplotlib

try:
   import cPickle as pickle
except:
   import pickle
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')


import UCSC_revised
import GenomeData
import SeparateByChrom
import associate_tags_with_regions
import get_total_tag_counts
import Utility_extended
import Entrez
from Entrez import EntrezGene
from Entrez import KnownEntrezGenes

plus = re.compile("\+");
minus = re.compile("\-");



def test_distribution_dic(dic, entrez_id):
	"""
	{entrezID:[((start, end), read_count)]}
	"""
	print "Entrez ID: ", entrez_id
	for element in dic[entrez_id]:
		print element


def calculateExonIntrons(entrez_genes, bedfile, column_index, chroms,  fragment_size, totalcount, out_file=None):
	"""
	entrez_genes is a EntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	
	return:
	all_reads_on_shared_exons = {} # {entrezID:[((start, end), read_count)]}
	all_reads_on_shared_introns = {} # {entrezID:[((start, end), read_count)]}
	all_reads_on_merged_transcripts = {} #{entrezID:[((start, end), read_count)]}
	all_summary = {} # {entrezID:{attribute:value}}
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
	lib_name = (bedfile).split('/')[-1] # remove directory
	suffix = lib_name.split('.')[-1] # txt
	lib_name = lib_name.split('.')[0] 
	extension = "-" + lib_name +'.' + suffix +"1"
	if Utility_extended.fileExists(bedfile):
		if Utility_extended.chrom_files_exist(chroms, extension) != 1:
			# Separate by chrom and sort by start
			print chroms, extension, " files do not exist, separate by chroms. "
			Utility_extended.separate_by_chrom_sort(chroms, bedfile, extension, [column_index])
	else:
		print bedfile, " is not found";
		sys.exit(1)
	
	all_reads_on_shared_exons = {} # {entrezID:[((start, end), read_count)]}
	all_reads_on_shared_introns = {} # {entrezID:[((start, end), read_count)]}
	all_reads_on_merged_transcripts = {} #{entrezID:[((start, end), read_count)]}
	all_summary = {} # {entrezID:{attributes}}
	
	for chrom in chroms:
		chrombed = chrom + extension
		if chrom in entrez_genes.chroms:
			entrez_genes_by_chrom =  Entrez.KnownEntrezGenes([chrom], entrez_genes.subset_by_chrom(chrom))
			(reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary) =  calculateExonIntrons_by_chrom (entrez_genes_by_chrom, chrombed, fragment_size, totalcount, out_file)
			#if chrom == chroms[0]:
				#myid = reads_on_shared_exons.keys()[0]
				#test(entrez_genes_by_chrom, reads_on_shared_introns, myid)
			all_reads_on_shared_exons.update(reads_on_shared_exons)
			all_reads_on_shared_introns.update(reads_on_shared_introns)
			all_reads_on_merged_transcripts.update(reads_on_merged_transcripts)
			all_summary.update(summary)
			print len(all_summary.keys())
		
	SeparateByChrom.cleanup(chroms, extension)
	return (all_reads_on_shared_exons, all_reads_on_shared_introns, all_reads_on_merged_transcripts, all_summary)

def calculateExonIntrons_by_chrom (entrez_genes_by_chrom, chrombed, fragment_size, totalcount, out_file=None):
	"""
	entrez_genes_by_chrom is a Entrez class object in a particular chrom
	bedfile is for a particular chrom
	totalcount: for RPKM
	out_file default is None and not writing to file
	return:
		reads_on_shared_exons={}	# {entrezID:[((start, end), read_count)]}
		reads_on_shared_introns={}	# {entrezID:[((start, end), read_count)]}
		reads_on_merged_transcripts={}	# {entrezID:[((start, end), read_count)]}
		summary = {} #	# {entrezID:{attribute:value}}
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
	if out_file is not None:
		outf = open(out_file, 'a')
	
	if Utility_extended.fileExists(chrombed):		
		# load in each read
		tag_position_list = []
		inf = open(chrombed,'r')
		for line in inf:
			if not re.match("#", line):
				line = line.strip()
				sline = line.split()
				tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
		inf.close()	
		if not Utility_extended.is_list_sorted(tag_position_list):
			tag_position_list.sort() #[tag_positions]
		
		tag_positions = [(item, 0) for item in tag_position_list] # convert its form to acceptable by get_read_counts_on_regions
		reads_on_shared_exons={}	# {entrezID:[((start, end), read_count)]}
		reads_on_shared_introns={}	# {entrezID:[((start, end), read_count)]}
		reads_on_merged_transcripts={}	# {entrezID:[((start, end), read_count)]}
		summary = {} # {entrezID:{attribute:value}}
			
		if entrez_genes_by_chrom.num_genes > 0:
			for entrez_id in entrez_genes_by_chrom.entrez_ids:
				
				gene = (entrez_genes_by_chrom.entrez_genes)[entrez_id]
				
				shared_exons = gene.shared_exonic_regions #sorted in absolute coordinate
				if shared_exons == []:# No shared extrons
					reads_on_shared_exons[entrez_id] = []
					shared_exons_rc = 0
					shared_exons_total_length = 0
					shared_exon_RPKM = 0
				else:
					# tag_positions is required to be [(position, annotation)], here annotation is not present
					# element in result_list has structure: [((start, end), read_count)]
					result_list = Utility_extended.get_read_counts_on_regions (tag_positions, shared_exons) 
					reads_on_shared_exons[entrez_id] = result_list
					
					shared_exons_total_length = sum([region[1] - region[0] + 1 for region in  shared_exons])
					shared_exons_rc = sum([item[1] for item in result_list])
					shared_exon_RPKM = shared_exons_rc * (1000.0/shared_exons_total_length) * (1000000.0/float(totalcount))
					
				shared_introns = gene.shared_intronic_regions #sorted
				if shared_introns == []: # No shared introns
					reads_on_shared_introns[entrez_id] = []
					shared_introns_rc = 0
					shared_introns_total_length = 0
					shared_intron_RPKM = 0
				else:
					# element in result_list has structure: [((start, end), read_count)]
					result_list = Utility_extended.get_read_counts_on_regions (tag_positions, shared_introns) 
					reads_on_shared_introns[entrez_id] = result_list
					shared_introns_total_length = sum([region[1]-region[0]+1 for region in  shared_introns])
					shared_introns_rc = sum([item[1] for item in result_list])
					shared_intron_RPKM = shared_introns_rc * (1000.0/shared_introns_total_length) * (1000000/float(totalcount))
					
				merged_exons = gene.merged_exonic_regions #sorted
				if merged_exons == []:
					merged_exons_total_length = 0
					merged_exons_rc = 0
					merged_exon_RPKM = 0.0
				else:
					result_list = Utility_extended.get_read_counts_on_regions(tag_positions, merged_exons)
					merged_exons_total_length = sum([region[1]-region[0]+1 for region in  merged_exons])
					merged_exons_rc = sum([item[1] for item in result_list])
					merged_exon_RPKM = merged_exons_rc * (1000.0/merged_exons_total_length) * (1000000/float(totalcount))
				
				merged_transcript = gene.boundaries #[(start, end)]
				# element in result_list has structure: [((start, end), read_count)]
				result_list = Utility_extended.get_read_counts_on_regions (tag_positions, merged_transcript)
				reads_on_merged_transcripts[entrez_id] = result_list
				merged_transcript_length = sum([region[1]-region[0]+1 for region in  merged_transcript])
				merged_transcript_rc = sum([item[1] for item in result_list])
				merged_transcript_RPKM = merged_transcript_rc * (1000.0/merged_transcript_length) * (1000000/float(totalcount))

				
				if out_file is not None:
					gene_symbol = []
					for transcript in gene.transcripts:
						if transcript.additional_annotations[0] not in gene_symbol:
							gene_symbol.append(transcript.additional_annotations[0])
					outline = str(entrez_id) + '\t' + str(merged_exons_rc) + '\t' + str(merged_exons_total_length) + '\t' + str(merged_exon_RPKM) + '\t' + str(shared_exons_rc) + '\t' + str(shared_exons_total_length) + '\t' + str(shared_exon_RPKM) + '\t' + str(shared_introns_rc) + '\t' + str(shared_introns_total_length) + '\t' +str(shared_intron_RPKM) + '\t' + str(merged_transcript_rc) + '\t' + str(merged_transcript_length) + '\t' + str(merged_transcript_RPKM) + '\t' + ','.join([transcript.name for transcript in gene.transcripts]) + '\t' + ','.join(gene_symbol) + '\n'
					outf.write(outline)
				
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
	
	if out_file is not None:
		outf.close()
	
	return (reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary)	
	
	
######################################################################################
######################################################################################
# intron level instead of gene level analysis of iri
######################################################################################
######################################################################################

def readin(mypickle):
	distribution = open(mypickle, 'rb')
	rc_distribution = pickle.load(distribution) 
	distribution.close()
	#print "The keys of " + opt.resting_intron_pickle + " are: ", rc_distribution.keys()
	print "There are ", len(rc_distribution.keys()), " genes in ", mypickle
	return rc_distribution

def test(reads_on_regions, my_entrez_id): 
	"""
	reads_on_regions: {entrezID:[((start, end), read_count)]}
	"""
	for item in (reads_on_regions)[my_entrez_id]:
		print item


def rpkms_on_regions(reads_on_regions, totalcount, id_set=None):
	"""
	Calculate the rpkms. 
	reads_on_regions:{entrezID:[((start, end), read_count)]}
	Return {entrezID:[((start, end), density)]}

	"""
	myids = Utility_extended.get_subset_ids_from_dic(reads_on_regions, id_set)

	rpkms={}
	for myid in myids:
		densities_on_gene = []
		for item in reads_on_regions[myid]:
			start = (item[0])[0]
			end = (item[0])[1]
			rc = item[1]
			density = float(rc) * (1000.0/(end-start+1)) * (1000000/float(totalcount))
			densities_on_gene.append((item[0], density))
		rpkms[myid] = densities_on_gene
	return rpkms

	
def get_expression_rpkms_on_shared_exons(rc_on_shared_exons, totalcount, id_set=None):
	"""
	rc_on_shared_exons:{entrezID:[((start, end), read_count)]}
	Returns {id:value}
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rc_on_shared_exons, id_set)
	
	expression_rpkm = {}
	for myid in myids:
		exons = rc_on_shared_exons[myid]
		
		shared_exons_total_length = sum([item[0][1]-item[0][0]+1 for item in  exons])
		shared_exons_rc = sum([item[1] for item in exons])
		expression_rpkm [myid] = shared_exons_rc * (1000.0/shared_exons_total_length) * (1000000/float(totalcount))
	return expression_rpkm


def select_by_expression_threshold(rc_on_shared_exons, totalcount, expression_cutoff, id_set=None):
	"""
	expression cutoff is in terms of rpkm.
	"""
	expression_rpkm = get_expression_rpkms_on_shared_exons(rc_on_shared_exons, totalcount, id_set)
	selected_ids=[]
	for myid in expression_rpkm.keys():
		expression = expression_rpkm[myid]
		if expression>=expression_cutoff:
			selected_ids.append(myid)
	return selected_ids

def get_trons_rc(rc_on_trons):
	"""
	get the total number of reads on (shared) introns or exons of a gene
	input: [((start, end), rc)]
	"""
	rc = sum([item[1] for item in rc_on_trons])
	return rc
	
def get_trons_length(rc_on_trons):
	"""
	get the total length of (shared) introns or exons of a gene
	input: [((start, end), rc)]
	"""
	rc = sum([(item[0][1]-item[0][0] + 1) for item in rc_on_trons])
	return rc
		
	
def get_iri_by_gene(rc_on_shared_introns, rc_on_shared_exons, rc_threshold, id_set=None):
	"""
	normalize the intron read density by the average exon read density of that gene
	
	rc_on_shared_introns: {entrezID:[((start, end), rc)]}
	rc_on_shared_exons: {entrezID:[((start, end), rc)]}
	
	expression cutoff is implemented by id_set and rc_threshold
	returns {entrezID:[((start, end), iri)]}
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rc_on_shared_introns, id_set)
	myids = Utility_extended.get_subset_ids_from_dic(rc_on_shared_exons, myids)
	
	iri_by_gene = {} #{entrez_id:[((start, end), iri)]}
	for entrez_id in myids:
		exons_rc = get_trons_rc(rc_on_shared_exons[entrez_id])
		if exons_rc > rc_threshold:
			iri_by_gene[entrez_id] = []
			exons_length = get_trons_length(rc_on_shared_exons[entrez_id])
			rc_density_on_exons = float(exons_rc)/exons_length
			for intron in rc_on_shared_introns[entrez_id]:
				intron_rc = intron[1]
				intron_coordinate =  intron[0]
				intron_length = intron_coordinate[1] - intron_coordinate[0] + 1
				rc_density_on_intron = float(intron_rc)/intron_length
				iri = rc_density_on_intron/rc_density_on_exons
				iri_by_gene[entrez_id].append((intron_coordinate, iri))
	
	return iri_by_gene	
	

def rank_iri_by_intron_length(iri, id_set=None):
	"""
	The module is to test whether long introns likely will have higher intron retention
	It uses the gene-based normalization for 
	
	iri: {entrezID:[((start, end), iri)]}
	returns [(id, length, iri_value)]
	"""
	#print "3"
	myids = Utility_extended.get_subset_ids_from_dic(iri, id_set)
	#print len(myids)
	
	ranked_list = []
	for myid in myids:
		#print "The id is ", myid
		introns = iri[myid]
		for item in introns:
			start = item[0][0]
			end = item[0][1]
			length = end - start + 1
			iri_value = item[1]
			ranked_list.append((myid, length, iri_value))
	ranked_list.sort(key = itemgetter(1))
	
	return ranked_list	
	
def log_transform_rpkm_distribution(rd_distribution, pc=0.0000000001, id_set=None):
	"""
	returns {entrezID:[((start, end), density)]
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rd_distribution, id_set)
				
	log_read_density_distribution = {}
	for entrez_id in myids:
		rpkms = rd_distribution[entrez_id]
		log_read_density_distribution[entrez_id] = [(item[0], log(item[1]+pc,2)) for item in rpkms]
	return log_read_density_distribution
	


def adjust_order_by_strand(rc_distribution, known_entrez_genes):
	"""
	Reorder the sequence of iri according to strand information so that one can aggregate by intron index
	return: {entrez_id: [((start, end), rc)]}
	"""
	my_distribution = {}
	myids = list(set(rc_distribution.keys()) & set(known_entrez_genes.entrez_ids))
	for entrez_id in myids:
		rcs = (rc_distribution)[entrez_id]
		gene = known_entrez_genes.entrez_genes[entrez_id]
		if minus.match(gene.strand):
			rcs.reverse()
		my_distribution[entrez_id] = copy.deepcopy(rcs)
	return my_distribution
				
				
def get_tron_rpkm_histogram(rpkm_distribution, mylabel, mytitle, id_set=None):
	"""
	A significant fraction (35%) of intronic regions are free of reads. Is it because of mappability? The simplest approach is the output these regions and check them out on genome browser. Another approach is to calculate an irrelevant ChIP-Seq (K36me3)/RNA-Seq library and compare the same regions, and use the scatter plot to find out the mappability.
	
	rpkm_distribution
	
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rpkm_distribution, id_set)
	myids.sort()# enable predictable behavior of entrez_id order
	
	all_list = [] # One entry each tronic regions
	
	for entrez_id in myids:
		rpkms = [item[1] for item in (rpkm_distribution)[entrez_id]]
		all_list.extend(rpkms)
		
	plt.clf()
	plt.hist(all_list, bins=40, color='r', normed=True, log=True )
	# mytitle = "Intron read density (normalized by average exon read density of respective gene) histogram"
	plt.title(mytitle)
	plt.xlabel(mylabel)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")		
	return 	all_list	

def get_fraction_retained_intron_per_gene_histogram(rc_distribution, rc_threshold, mytitle, id_set=None):
	"""
	rc_distribution:{entrezID:[((start, end), rc)]
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rc_distribution, id_set)
	myids.sort() # enable predictable behavior of chrom order
	
	all_list = {} # one entry each gene
	total = 0.0
	total_above_threshold = 0.0 
	total_no_tron_genes = 0
	for entrez_id in myids:
		rcs = [item[1] for item in (rc_distribution)[entrez_id]]#read counts[]
		if len(rcs) > 0:
			above_threshold =0
			for item in rcs:
				if item >= rc_threshold:
					above_threshold += 1
			above_threshold_fraction = float(above_threshold)/(len(rcs)) 
			all_list[entrez_id] = above_threshold_fraction
			total +=len(rcs)
			total_above_threshold += above_threshold
		else:
			#print "%s has no trons" %entrez_id
			total_no_tron_genes += 1
	print "%d out of %d genes have no trons" %(total_no_tron_genes	, len(myids))	
	print "The number of introns with read count >= %d is %d, fraction  %f" %(rc_threshold,total_above_threshold, float(total_above_threshold)/total)
	plt.clf()
	plt.hist(all_list.values(), bins=100, color='r', normed=True )
	# mytitle = "Intron read density (normalized by average exon read density of respective gene) histogram"
	plt.title(mytitle)
	plt.xlabel("Fraction of retained introns per gene ")
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")		
	return all_list		

	
def get_rd_fluctuation_histogram(rd_distribution, num_introns_cutoff, mylabel, mytitle, pc = 0.000000001, id_set=None):
	"""
	rd_distribution:{entrezID:[((start, end), rd)]}
	For each simple gene, calculate fluctuation of per-intron read density normalized by expression: standard deviation {entrez_id: rf}
	
	can also be used for iri
	Return {entrez_id:value}
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rd_distribution, id_set)
	
	fluctuation = {}
	for entrez_id in myids:
			this_gene = rd_distribution[entrez_id]
			if len(this_gene) >= num_introns_cutoff:
				iris = [item[1] for item in this_gene]
				std = numpy.std(iris)
				fluctuation[entrez_id] = std
	plt.clf()
	plt.figure(1)
	plt.subplot(211)
	plt.hist(fluctuation.values(), bins=50, color='g', normed=True, log=True)
	plt.title(mytitle)
	plt.xlabel(mylabel)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	
	plt.subplot(212)
	plt.hist([log(item + pc, 2) for item in fluctuation.values()], bins=50, color='r', normed=True, log=True)
	plt.xlabel("log "+ mylabel)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")
	
	return fluctuation
	

def get_rd_relative_fluctuation_histogram(rd_distribution, num_introns_cutoff, mylabel, mytitle, pc = 0.000000001, id_set=None):
	"""
	rd_distribution:{entrezID:[((start, end), rd)]}
	For each simple gene, calculate a relative fluctuation of per-intron iri: standard deviation/mean {entrez_id: rf}
	The reason for relative fluctuation is for 
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rd_distribution, id_set)
	
	relative_fluctuation = {}
	for entrez_id in myids:
			this_gene = rd_distribution[entrez_id]
			if len(this_gene) >= num_introns_cutoff:
				iris = [item[1] for item in this_gene]
				mean = numpy.average(iris)
				if mean > 0:
					std = numpy.std(iris)
					rf = std/mean
					relative_fluctuation[entrez_id] = rf
	plt.clf()
	plt.figure(1)
	plt.subplot(211)
	plt.hist(relative_fluctuation.values(), bins=50, color='g', normed=True, log=True)
	plt.title(mytitle)
	plt.xlabel(mylabel)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	
	plt.subplot(212)
	plt.hist([log(item + pc, 2) for item in relative_fluctuation.values()], bins=50, color='r', normed=True, log=True)
	plt.xlabel("log "+ mylabel)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")
	
	return relative_fluctuation
	
	
######################################################################################
# intron level instead of gene level analysis of iri, only for simple genes
######################################################################################	
	
def get_intron_level_iri_distribution(reads_on_shared_introns, reads_on_shared_exons, id_set, rc_threshold):
	"""
	Calculate IRI by the neighbouring exons, instead of averaging over all the exons of the genes
	
	Simple genes are genes whose entrez id has only on refseq id. Only for these genes, we are sure that the introns are in between exons. Simple genes are reflected in id_set
	
	rc_threshold: if rc in its neighbouring exons is below cutoff, don't count
	
	returns {entrezID:[((start, end), iri)]}
	"""
	
	iri = {}

	mykeys = list( set((reads_on_shared_introns).keys()) & set(id_set))
	for entrez_id in mykeys:
		intron_reads = (reads_on_shared_introns)[entrez_id]
		exon_reads = (reads_on_shared_exons)[entrez_id]
		assert (len(exon_reads) == len(intron_reads) + 1)
		(iri)[entrez_id] = []
		for i in xrange(len(intron_reads)):
			coordinate = intron_reads[i][0]
			intron_length = coordinate[1] - coordinate[0] + 1
			intron_rc = intron_reads[i][1]
			intron_density = intron_rc/float(intron_length)
			exon1_length = exon_reads[i][0][1] - exon_reads[i][0][0] + 1
			exon2_length = exon_reads[i+1][0][1] - exon_reads[i+1][0][0] + 1
			neighbor_exons_rc = exon_reads[i][1] + exon_reads[i+1][1]
			neighbor_exons_density = float(neighbor_exons_rc)/(exon1_length + exon2_length)
			if neighbor_exons_rc > rc_threshold:
				if intron_rc == 0:
					myiri = 0
				else:
					myiri = intron_density/neighbor_exons_density
				(iri)[entrez_id].append((coordinate, myiri))
	return iri

def get_distribution_ito_tron_number(rc_distribution, id_set=None):
	"""
	Most appropriate for simple genes
	ito: in terms of 
	Resection rc_distribution so that it is organized by the number of trons 
	Return: {number_of_trons:{entrez_id:[((start, end), iri)]}}
	"""
	myids = Utility_extended.get_subset_ids_from_dic(rc_distribution, id_set)

	distribution_ito_tron_number = {} # {number_of_trons:{entrez_id:[((start, end), iri)]}}
	for entrez_id in myids:
		rcs = (rc_distribution)[entrez_id]
		number_of_trons = len(rcs)
		if number_of_trons not in distribution_ito_tron_number.keys():
			distribution_ito_tron_number[number_of_trons] = {}
		else:
			(distribution_ito_tron_number[number_of_trons])[entrez_id] = rcs
	return distribution_ito_tron_number	

def Internal_vs_external_trons_boxplot(distribution_ito_tron_number, name, threshold=3, pc = 0.0001):
	"""
	Most appropriate for simple genes
	threshold: number of intron threshold
	distribution_ito_tron_number comes from above: {number_of_trons:{entrez_id:[((start, end), iri)]}}
	Return [[first intron iris],[middle intron iris],[last intron iris]]
	"""
	my_list = []
	for i in xrange(3):
		my_list.append([])
	number_of_trons = sorted(distribution_ito_tron_number.keys())
	for i in number_of_trons:
		if i >= threshold:
			for entrez_id in distribution_ito_tron_number[i].keys():
				# First intron
				(my_list[0]).extend([((distribution_ito_tron_number[i])[entrez_id])[0][1]]) # the argument of extend is a list
				# middle introns
				temp = ((distribution_ito_tron_number[i])[entrez_id])[1:-1] #[((start, end), iri)]
				(my_list[1]).extend([item[1] for item in temp])
				# last intron
				(my_list[2]).extend([((distribution_ito_tron_number[i])[entrez_id])[-1][1]])
	plt.clf()
	plt.boxplot(my_list)
	plt.savefig(name + ".png", format="png")
	return my_list	

def generate_characteristic_figures(distribution_ito_tron_number, m, name, pc = 0.001):
	"""
	distribution_ito_tron_number: {number_of_introns:{id:[((start,end), rc)]}}
	m: the number of introns
	name: output file name
	Most appropriate for simple genes
	Draw heat map for a list of genes with equal number of introns.
	
	distribution_ito_tron_number:get_distribution_ito_tron_number
	"""
	my_list = []
	if ( m in distribution_ito_tron_number.keys()):
		# Get the relevant rpkms
		for entrez_id in (distribution_ito_tron_number[m]).keys():
			rcs = [item[1] for item in (distribution_ito_tron_number[m])[entrez_id]]
			my_list.append([log(item+pc, 2) for item in rcs])
		print "There are ", len(my_list), " genes that has ", m, " read counts"  
	
		my_array = numpy.array(my_list)
		my_transposed_array = numpy.transpose(my_array)
		
		# draw boxplot for sequential distribution of the columns.
		plt.clf()
		plt.boxplot(my_array)
		plt.savefig(name + "_" + str(m) + "_boxplot.png", format="png")
	
		# heatmap

		r.library("gplots")
		
		my_special_list = []
		for i in xrange(len(my_list)):
			for item in my_list[i]:
				my_special_list.append(float(item))
		numeric=r['is.numeric']
		asnumeric=r['as.numeric']		
		heatmap2 = r['heatmap.2']
		
		#print numeric(asnumeric(my_special_list))
		
		matrix = robjects.r['matrix'](asnumeric(my_special_list), ncol=m, nrow=len(my_list))	
		#matrix = r.matrix(rlist, ncol=m, nrow=len(my_list))	
		#print r.dim(matrix)
		
		grdevices = importr('grDevices')
		grdevices.png(file = name + "_" + str(m) + "_clustered-R.png", width=512, height=512)
		#heatmap2(matrix, Colv=False, Rowv=True, dendrogram="row", key=True, keysize=1, col=r.redgreen(75), trace="none", scale="row", main="Euclidean")	
		heatmap2(matrix, Colv=False, Rowv=True, dendrogram="row", key=True, keysize=1, col=r.redgreen(75), trace="none", scale="row", main="Euclidean")	
		grdevices.dev_off()
	
	return distribution_ito_tron_number

	
	
def cell_type_sepecific_analysis(reads_on_shared_introns, reads_on_shared_exons, entrez_genes, totalcount, cell_type,  expression_cutoff, pc):
	"""
	cell_type: "resting" or "activated"
	"""
	# enforce the expression threshold
	rpkms_on_shared_exons = rpkms_on_regions(reads_on_shared_exons, totalcount)
	selected_ids = select_by_expression_threshold(rpkms_on_shared_exons, totalcount, expression_cutoff)

	# Calculate the rpkms_on_shared_introns
	rpkms_on_shared_introns = rpkms_on_regions(reads_on_shared_introns, totalcount, selected_ids)
	# Normalize the intron_read_density by the exons-rpkm of the respective gene
	rc_threshold = 10
	iri_by_gene = get_iri_by_gene(reads_on_shared_introns, reads_on_shared_exons, rc_threshold, selected_ids) #{entrezID:[((start, end), iri)]}
	
	mylabel = "per gene relative fluctuation of intron read density"
	mytitle = "Distribution_of_per_gene_relative_fluctuation_of_intron_read_density_" + cell_type
	num_introns_cutoff = 4
	read_density_relative_fluctuation = get_rd_fluctuation_histogram (iri_by_gene, num_introns_cutoff, mylabel, mytitle, pc) #{entrez_id:value}
	log_iri_by_gene = log_transform_rpkm_distribution(iri_by_gene, pc) #{entrezID:[((start, end), density)]
	mylabel = "log(Normalized read density,2)"
	mytitle = "Distribution_of_intronic_rpkms_normalized_by_expression"+ cell_type
	get_tron_rpkm_histogram(log_iri_by_gene, mylabel, mytitle)
	
	rc_threshold = 10
	mytitle = "Distribution_of_genes_in_terms_of_fraction_of_introns_with_rc_above_" + str(rc_threshold) + cell_type
	get_fraction_retained_intron_per_gene_histogram(reads_on_shared_introns, rc_threshold, mytitle)
	
	my_intron_dic = get_distribution_ito_tron_number(iri_by_gene)
	my_intron_log_dic = get_distribution_ito_tron_number(log_iri_by_gene)
	name = "Internal_vs_external_introns_boxplot_" + cell_type
	#[[first intron iris], [middle intron iris], [last intron iris]]
	internal_vs_external_introns_iris = Internal_vs_external_trons_boxplot(my_intron_log_dic, name, 4)
	#for number_of_trons in xrange(4,10):
		#generate_characteristic_figures(my_intron_dic, number_of_trons, "intron_read_density_" + str(number_of_trons) + "_" + opt.name)

	return (read_density_relative_fluctuation, internal_vs_external_introns_iris)
			
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--forwardreadfile", action="store", type="string", dest="ReadsOnForwardStrand", help="input bed file for RNASeq raw reads on forward strand", metavar="<file>")
	parser.add_option("-r", "--reversereadfile", action="store", type="string", dest="ReadsOnReverseStrand", help="input bed file for RNASeq raw reads on reverse strand", metavar="<file>")
	parser.add_option("-g", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determines the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-u", "--entrez_genes_file", action="store", type="string", dest="entrez_genes", metavar="<file>", help="file with curated known genes clustered by entrez ID in pickle format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name for genes and tag numbers")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	test=0
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	
	startTime = time.time()
	
	##################################################################3
	#The column numbers are 1 based instead of 0 based!
	#For positive strand
	start_index_P = 2
	#For negative strand
	start_index_N = 3
	##################################################################3
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	# entrez_gene_collection is a EntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
	
		
	# test module
	test = 0
	if test == 1:
		print "Testing gene structure"
		test_id = 54
		Entrez.test_gene_structure(entrez_gene_collection, test_id)

	
	totalcount_F = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnForwardStrand);
	totalcount_R = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnReverseStrand);
	totalcount = totalcount_F + totalcount_R
	print totalcount_F, totalcount_R

	#Clear the file.
	outf = open(opt.outfile, 'w')
	outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	outf.write(outline)
	outf.close()
	
	# The RNA seq data are strand specific. Only use + reads on genes on forward strand, and - reads on genes on reverse strand.
	print "Process genes on forward strand"
	entrez_ids_on_forward_strand = entrez_gene_collection.get_strand_specific_ids("+")
	print "There are ", len(entrez_ids_on_forward_strand), " Entrez IDs on forward strand."
	entrez_gene_subset = Entrez.KnownEntrezGenes(chroms, entrez_gene_collection.subset(entrez_ids_on_forward_strand))
	
	(forward_reads_on_shared_exons, forward_reads_on_shared_introns, forward_reads_on_merged_transcripts, forward_summary) = calculateExonIntrons(entrez_gene_subset, opt.ReadsOnForwardStrand, start_index_P, chroms, opt.fragment_size, totalcount, opt.outfile)
	
	print "Process genes on reverse strand"
	entrez_ids_on_reverse_strand = entrez_gene_collection.get_strand_specific_ids("-")
	print "There are ", len(entrez_ids_on_reverse_strand), " Entrez IDs on reverse strand."
	entrez_gene_subset = Entrez.KnownEntrezGenes(chroms, entrez_gene_collection.subset(entrez_ids_on_reverse_strand))
	
	(reverse_reads_on_shared_exons, reverse_reads_on_shared_introns, reverse_reads_on_merged_transcripts, reverse_summary) = calculateExonIntrons(entrez_gene_subset, opt.ReadsOnReverseStrand, start_index_N, chroms, opt.fragment_size, totalcount, opt.outfile)
	
	#combine the densities
	# {entrezID:[((start, end), read_count)]}
	reads_on_shared_exons = {}
	reads_on_shared_exons.update(forward_reads_on_shared_exons)
	reads_on_shared_exons.update(reverse_reads_on_shared_exons)
	name = opt.outfile + "_shared_exons.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_shared_exons, output)
	output.close()	
	
	if test == 1:
		test_distribution_dic(reads_on_shared_exons, test_id)
	
	
	# {entrezID:[((start, end), read_count)]}
	reads_on_shared_introns = {}
	reads_on_shared_introns.update(forward_reads_on_shared_introns)
	reads_on_shared_introns.update(reverse_reads_on_shared_introns)
	#store the info in a pickle file
	name = opt.outfile + "_shared_introns.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_shared_introns, output)
	output.close()
	
	if test == 1:
		test_distribution_dic(reads_on_shared_introns, test_id)
	
	reads_on_merged_transcripts = {}
	reads_on_merged_transcripts.update(forward_reads_on_merged_transcripts)
	reads_on_merged_transcripts.update(reverse_reads_on_merged_transcripts)
	#store the info in a pickle file
	name = opt.outfile + "_merged_transcripts.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_merged_transcripts, output)
	output.close()
	
	summary = {}
	summary.update(forward_summary)
	summary.update(reverse_summary)
	name = opt.outfile + "_summary.pkl"
	output = open(name, 'wb')
	pickle.dump(summary, output)
	output.close()
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)