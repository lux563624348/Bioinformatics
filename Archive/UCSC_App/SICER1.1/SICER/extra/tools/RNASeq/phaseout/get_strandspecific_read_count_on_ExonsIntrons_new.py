#!/usr/bin/env python

# Calculate strand-specific read counts on intronic, exonic regions with strand specific reads (eg, strand-specific RNA Seq)

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
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/JunZhu/modules')

import UCSC_revised
import GenomeData
import SeparateByChrom
import associate_tags_with_regions
import get_total_tag_counts
import Utility_extended
import Entrez
from Entrez import EntrezGene
from Entrez import KnownEntrezGenes
import AnalyzeRNASeqPASeqChIPSeq

plus = re.compile("\+");
minus = re.compile("\-");



def test_distribution_dic(dic, entrez_id):
	"""
	{entrezID:[((start, end), read_count)]}
	"""
	print "Entrez ID: ", entrez_id
	for element in dic[entrez_id]:
		print element


def calculateExonIntrons(entrez_genes, bedfile, column_index, chroms,  fragment_size, totalcount, out_file):
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
	all_summary = {}
	
	for chrom in chroms:
		chrombed = chrom + extension
		entrez_genes_by_chrom =  Entrez.KnownEntrezGenes([chrom], entrez_genes.subset_by_chrom(chrom))
		(reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary) =  calculateExonIntrons_by_chrom (entrez_genes_by_chrom, chrombed, fragment_size, totalcount, out_file)
		#if chrom == chroms[0]:
			#myid = reads_on_shared_exons.keys()[0]
			#test(entrez_genes_by_chrom, reads_on_shared_introns, myid)
		all_reads_on_shared_exons.update(reads_on_shared_exons)
		all_reads_on_shared_introns.update(reads_on_shared_introns)
		all_reads_on_merged_transcripts.update(reads_on_merged_transcripts)
		all_summary.update(summary)
		
	SeparateByChrom.cleanup(chroms, extension)
	return (all_reads_on_shared_exons, all_reads_on_shared_introns, all_reads_on_merged_transcripts, summary)

def calculateExonIntrons_by_chrom (entrez_genes_by_chrom, chrombed, fragment_size, totalcount, out_file):
	"""
	entrez_gene_collection is a class object in a particular chrom
	bedfile is for a particular chrom

	"""
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
		summary = {} #	# {entrezID:{attribute_name:value}}
			
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

				gene_symbol = []
				for transcript in gene.transcripts:
					if transcript.additional_annotations[0] not in gene_symbol:
						gene_symbol.append(transcript.additional_annotations[0])
				outline = str(entrez_id) + '\t' + str(merged_exons_rc) + '\t' + str(merged_exons_total_length) + '\t' + str(merged_exon_RPKM) + '\t' + str(shared_exons_rc) + '\t' + str(shared_exons_total_length) + '\t' + str(shared_exon_RPKM) + '\t' + str(shared_introns_rc) + '\t' + str(shared_introns_total_length) + '\t' +str(shared_intron_RPKM) + '\t' + str(merged_transcript_rc) + '\t' + str(merged_transcript_length) + '\t' + str(merged_transcript_RPKM) + '\t' + ','.join([transcript.name for transcript in gene.transcripts]) + '\t' + ','.join(gene_symbol) + '\n'
				outf.write(outline)
				
				summary[entrez_id] = {}
				(summary[entrez_id])["merged_exons_rc"] = merged_exons_rc
				(summary[entrez_id])["merged_exons_RPKM"] = merged_exon_RPKM
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
	outf.close()
	
	return (reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary)
	

######################################################################################
# intron level (instead of gene level) analysis of iri
######################################################################################
def readin(mypickle):
	distribution = open(mypickle, 'rb')
	rc_distribution = pickle.load(distribution) 
	distribution.close()
	#print "The keys of " + opt.resting_intron_pickle + " are: ", rc_distribution.keys()
	print "There are ", rc_distribution.keys(), " genes in ", mypickle
	return rc_distribution

def test(reads_on_regions, my_entrez_id): 
	"""
	reads_on_regions: {entrezID:[((start, end), read_count)]}
	"""
	for item in (reads_on_regions)[entrez_id]:
		print item

def rpkms_on_regions(reads_on_regions, totalcount, id_set=None):
	"""
	Calculate the rpkms. Return {entrezID:[((start, end), density)]}
	Still different from rpkm
	"""
	if id_set is None:
		myids = reads_on_regions.keys()
	else:
		myids = list(set(reads_on_regions.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	rpkms={}
	for myid in myids:
		densities_on_gene = []
		for item in reads_on_regions[myid]:
			start = (item[0])[0]
			end = (item[0])[1]
			rc = item[1]
			density = float(rc) * (1000.0/(end-start+1)) * (1000000/float(totalcount))
			densities_on_gene.append(item[0], density)
		rpkms[myid] = densities_on_gene
	return rpkms

def get_expression_rpkms_on_shared_exons(rc_on_shared_exons, totalcount, id_set=None):
	"""
	Returns {id:value}
	"""
	if id_set is None:
		myids = rc_on_shared_exons.keys()
	else:
		myids = list(set(rc_on_shared_exons.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
	expression_rpkm = {}
	for myid in rc_on_shared_exons:
		exons = rc_on_shared_exons[myid]
		
		shared_introns_total_length = sum([region[1]-region[0]+1 for region in  exons])
		shared_introns_rc = sum([item[1] for item in exons])
		expression_rpkm [myid] = shared_introns_rc * (1000.0/shared_introns_total_length) * (1000000/float(totalcount))
	return expression_rpkm


def select_by_expression_threshold(rc_on_shared_exons, totalcount, expression_cutoff, id_set=None):
	"""
	expression cutoff is in terms of rpkm.
	"""
	if id_set is None:
		myids = rpkms_on_shared_exons.keys()
	else:
		myids = list(set(rpkms_on_shared_exons.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
	expression_rpkm = get_expression_rpkms_on_shared_introns(rc_on_shared_exons, totalcount, id_set)
	selected_ids=[]
	for myid in expression_rpkm.keys():
		expression = expression_rpkm[myid]
		if expression>=expression_cutoff:
			selected_ids.append(myid)
	return selected_ids

def normalize_intron_read_density(rc_on_shared_introns, rc_on_shared_exons, totalcount, id_set=None):
	"""
	normalize the intron read density by the average exon read density of that gene, which is available from summary
	expression cutoff is implemented by id_set
	returns {entrezID:[((start, end), density)]}
	"""
	if id_set is None:
		myids = list(set(rc_on_shared_introns.keys()) & set(rc_on_shared_introns.keys()))
	else:
		myids = list(set(rpkms_on_shared_introns.keys()) & set(rc_on_shared_introns.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
	rpkm_on_shared_introns = rpkms_on_regions(rc_on_shared_introns, totalcount, myids)
	expression_rpkm = get_expression_rpkms_on_shared_exons(rc_on_shared_exons, totalcount, myids)
	
	normalized_intron_read_density = {} #{entrez_id:[((start, end), rpkm)]}
	for entrez_id in myids:
		intron_rpkms = rpkms_on_shared_introns[entrez_id]
		normalization = float(expression_rpkm[entrez_id])
		normalized_intron_read_density[entrez_id] = [(item[0], item[1]/normalization) for item in intron_rpkms]
	return normalized_intron_read_density


def rank_normalized_intron_read_density_by_intron_length(normalized_intron_read_density, id_set=None):
	"""
	The module is to test whether long introns likely will have higher intron retention
	It uses the gene-based normalization for 
	returns [(id, length, normalized_intron_read_density)]
	"""
	
	if id_set is None:
		myids = normalized_intron_read_density.keys()
	else:
		myids = list(set(normalized_intron_read_density.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))

	ranked_list = []
	for myid in myids:
		items = normalized_intron_read_density[myid]
		for item in items:
			start = item[0][0]
			end = item[0][1]
			length = end - start + 1
			normalized_intron_rd = item[1]
			ranked_list.append((myid, length, normalized_intron_rd))
	ranked_list.sort(key.itemgetter(1))
	
	#top vs bottom
	top = [item[2] for item in ranked_list[-1000:]]
	bottom = [item[2] for item in ranked_list[:1000]]
	my_list = [top, bottom]
	plt.clf()
	plt.boxplot(my_list)
	name = "iri_longest_1000_introns_vs_shortest_1000_introns_gene_based_normalization"
	plt.savefig(name + ".png", format="png")
	
	return ranked_list


def log_transform_rpkm_distribution(rd_distribution, pc=0.0000000001, id_set=None):
	"""
	returns {entrezID:[((start, end), density)]
	"""
	if id_set is None:
		myids = rd_distribution.keys()
	else:
		myids = list(set(rd_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
				print "%d ids are thrown out" % (len(id_set) - len(myids))
				
	log_read_density_distribution = {}
	for entrez_id in myids:
		rpkms = rd_distribution[entrez_id]
		log_read_density_distribution[entrez_id] = [(item[0], log(item[1]+pc,2)) for item in rpkms]
	return log_read_density_distribution

def adjust_order_by_strand(rc_distribution, known_entrez_genes):
	"""
	
	Reorder the sequence of iri according to strand information so that one can aggregate by intron index
	"""
	my_distribution = {}
	for entrez_id in (rc_distribution).keys():
		rcs = (rc_distribution)[entrez_id]
		gene = known_entrez_genes.entrez_genes[entrez_id]
		if minus.match(gene.strand):
			rcs.reverse()
		my_distribution[entrez_id] = rcs
	return my_distribution
				
				
def get_tron_rpkm_histogram(rpkm_distribution, mylabel, mytitle, id_set=None):
	"""
	A significant fraction (35%) of intronic regions are free of reads. Is it because of mappability? The simplest approach is the output these regions and check them out on genome browser. Another approach is to calculate an irrelevant ChIP-Seq (K36me3)/RNA-Seq library and compare the same regions, and use the scatter plot to find out the mappability.
	
	"""
	if id_set is None:
		myids = rpkm_distribution.keys()
	else:
		myids = list(set(rpkm_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
	all_list = [] # One entry each tronic regions
	myids.sort()# enable predictable behavior of entrez_id order
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

def get_fraction_retained_intron_per_gene_histogram(rc_distribution, rc_cutoff, mytitle, id_set=None):
	if id_set is None:
		myids = rc_distribution.keys()
	else:
		myids = list(set(rc_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	myids.sort() # enable predictable behavior of chrom order
	
	all_list = {} # one entry each gene
	total = 0.0
	total_non_zero = 0.0 
	for entrez_id in myids:
		rcs = [item[1] for item in (rc_distribution)[entrez_id]]
		non_zero =0
		for item in rcs:
			if item >= rc_cutoff:
				non_zero += 1
		non_zero_fraction = float(non_zero)/(len(rcs)) 
		all_list[entrez_id] = non_zero_fraction
		total +=len(rcs)
		total_non_zero += non_zero
	print "The overall fraction of introns with read count >= %d :  %d" %(rc_cutoff, total_non_zero/total)
	plt.clf()
	plt.hist(all_list.values(), bins=100, color='r', normed=True )
	# mytitle = "Intron read density (normalized by average exon read density of respective gene) histogram"
	plt.title(mytitle)
	plt.xlabel("Fraction of retained introns per gene ")
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")		
	return all_list		
	
def get_rd_fluctuation_histogram(rd_distribution, num_introns_cutoff, id_set=None):
	"""
	For each simple gene, calculate fluctuation of per-intron read density normalized by expression: standard deviation {entrez_id: rf}
	
	can also be used for iri
	"""
	if id_set is None:
		myids = rd_distribution.keys()
	else:
		myids = list(set(rd_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
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
	plt.xlabel(mylable)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	
	plt.subplot(212)
	plt.hist([log(item + pc, 2) for item in fluctuation.values()], bins=50, color='r', normed=True, log=True)
	plt.xlabel("log "+ mylable)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")
	
	
	return fluctuation

def get_rd_relative_fluctuation_histogram(rd_distribution, num_introns_cutoff, mylable, mytitle, pc = 0.000000001, id_set=None):
	"""
	For each simple gene, calculate a relative fluctuation of per-intron iri: standard deviation/mean {entrez_id: rf}
	The reason for relative fluctuation is for 
	"""
	if id_set is None:
		myids = rd_distribution.keys()
	else:
		myids = list(set(rd_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))
	
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
	plt.xlabel(mylable)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	
	plt.subplot(212)
	plt.hist([log(item + pc, 2) for item in relative_fluctuation.values()], bins=50, color='r', normed=True, log=True)
	plt.xlabel("log "+ mylable)
	plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	plt.savefig(mytitle + ".png", format="png")
	
	return relative_fluctuation
				

######################################################################################
# intron level instead of gene level analysis of iri, only for simple genes
######################################################################################

def find_simple_gene_iri_per_intron_distribution(reads_on_shared_introns, reads_on_shared_exons, id_set, rc_cutoff, pc = 0.5):
	"""
	Simple genes are genes whose entrez id has only on refseq id. Only for these genes, we are sure that the introns are in between exons. Simple genes are reflected in id_set
	
	rc_cutoff: if rc in an intron is below cutoff, don't count
	
	returns {entrezID:[((start, end), iri)}
	"""
	total = 0
	iri_distribution = {}

	mykeys = list( set((reads_on_shared_introns).keys()) & set(id_set))
	for entrez_id in mykeys:
		intron_reads = (reads_on_shared_introns)[entrez_id]
		exon_reads = (reads_on_shared_exons)[entrez_id]
		assert (len(exon_reads) == len(intron_reads) + 1)
		(iri_distribution)[entrez_id] = []
		for i in xrange(len(intron_reads)):
			coordinate = intron_reads[i][0]
			intron_rc = intron_reads[i][1]
			if intron_rc >= rc_cutoff:
				
				total += 1
				neighbor_exons_rc = (exon_reads[i][1] + exon_reads[i+1])/2.0
				if intron_rc == 0:
					myiri = 0
				else:
					myiri = (intron_rc + pc) / (neighbor_exons_rc + pc)
				(iri_distribution)[entrez_id].append((coordinate, myiri))
	return iri_distribution

def rank_iri_by_intron_length(iri_distribution, id_set=None):
	"""
	The module is to test whether long introns likely will have higher intron retention
	It uses the neighbor-exons normalization 
	returns [(id, length, iri)] sorted according to intron_length
	"""
	
	if id_set is None:
		myids = iri_distribution.keys()
	else:
		myids = list(set(iri_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))

	ranked_list = []
	for myid in myids:
		items = iri_distribution[myid]
		for item in items:
			start = item[0][0]
			end = item[0][1]
			length = end - start + 1
			iri = item[1]
			ranked_list.append((myid, length, iri))
	ranked_list.sort(key.itemgetter(1))
	
	#top vs bottom
	top = [item[2] for item in ranked_list[-1000:]]
	bottom = [item[2] for item in ranked_list[:1000]]
	my_list = [top, bottom]
	plt.clf()
	plt.boxplot(my_list)
	name = "iri_longest_1000_introns_vs_shortest_1000_introns_enighborexons_based_normalization"
	plt.savefig(name + ".png", format="png")
	return ranked_list

def get_distribution_ito_tron_number(rc_distribution, id_set=None):
	"""
	Most appropriate for simple genes
	
	ito: in terms of 
	Resection rc_distribution so that it is organized by the number of trons 
	"""
	if id_set is None:
		myids = rc_distribution.keys()
	else:
		myids = list(set(rc_distribution.keys()) & set(id_set))
		if len(myids) < len(id_set):
			print "%d ids are thrown out" % (len(id_set) - len(myids))

	distribution_ito_tron_number = {} # {number_of_trons:{entrez_id:[]}}
	for entrez_id in myids:
		rcs = (rc_distribution)[entrez_id]
		number_of_trons = len(rcs)
		if number_of_trons not in distribution_ito_tron_number.keys():
			distribution_ito_tron_number[number_of_trons] = {}
		else:
			(distribution_ito_tron_number[number_of_trons])[entrez_id] = rcs
	return distribution_ito_tron_number


def Internal_vs_external_trons_boxplot(distribution_ito_tron_number, name, num_intron_threshold=3, pc = 0.0001):
	"""
	Most appropriate for simple genes
	
	distribution_ito_tron_number comes from above: get_distribution_ito_tron_number
	
	"""
	my_list = []
	for i in xrange(3):
		my_list.append([])
	number_of_trons = sorted(distribution_ito_tron_number.keys())
	for i in number_of_trons:
		if i >= num_intron_threshold:
			for entrez_id in distribution_ito_tron_number[i].keys():
				# First intron
				(my_list[0]).extend([((distribution_ito_tron_number[i])[entrez_id])[0]]) # the argument of extend is a list
				# middle introns
				(my_list[1]).extend(((distribution_ito_tron_number[i])[entrez_id])[1:-1])
				# last intron
				(my_list[2]).extend([((distribution_ito_tron_number[i])[entrez_id])[-1]])
	#my_log_list = []
	#print len(my_list)
	#for i in xrange(3):
		#items = my_list[i]
		#print items[:10]
		#my_log_list.append ([log(item+pc, 2) for item in items])
	
	plt.clf()
	plt.boxplot(my_list)
	plt.savefig(name + ".png", format="png")
	return my_list

def generate_characteristic_figures(distribution_ito_tron_number, num_introns, name, pc = 0.001):
	"""
	Most appropriate for simple genes
	Draw heat map for a list of genes with equal number of introns.
	
	distribution_ito_tron_number:get_distribution_ito_tron_number
	"""
	my_list = []
	if ( num_introns in distribution_ito_tron_number.keys()):
		# Get the relevant rpkms
		for entrez_id in (distribution_ito_tron_number[num_introns]).keys():
			rcs = [item[1] for item in (distribution_ito_tron_number[num_introns])[entrez_id]]
			my_list.append([log(item+pc, 2) for item in rcs])
		print "There are ", len(my_list), " genes that has ", num_introns, " read counts"  
	
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


def cell_type_sepecific_analysis(reads_on_shared_introns, reads_on_shared_exons, entrez_genes, totalcount, cell_type expression_cutoff, pc):
	"""
	cell_type: "resting" or "activated"
	"""
	# enforce the expression threshold
	rpkms_on_shared_exons = rpkms_on_regions(reads_on_shared_exons, totalcount)
	selected_ids = select_by_expression_threshold(rpkms_on_shared_exons, expression_cutoff)
	
	# Calculate the rpkms_on_shared_introns
	rpkms_on_shared_introns = rpkms_on_regions(reads_on_shared_introns, totalcount, selected_ids)
	# Normalize the intron_read_density by the exons-rpkm of the respective gene
	normalized_rpkms_on_shared_introns = normalize_intron_read_density(rpkms_on_shared_introns, summary)
	
	
	mylable = "per gene relative fluctuation of intron read density"
	mytitle = "Distribution_of_per_gene_relative_fluctuation_of_intron_read_density_" + cell_type
	num_introns_cutoff = 4
	read_density_relative_fluctuation = get_rd_fluctuation_histogram (normalized_rpkms_on_shared_introns, num_introns_cutoff, mylable, mytitle)
	log_dist = log_transform_rpkm_distribution(normalized_rpkms_on_shared_introns, pc)
	mylabel = "log(Normalized read density,2)"
	mytitle = "Distribution_of_intronic_rpkms_normalized_by_expression"+ cell_type
	get_tron_rpkm_histogram(log_dist, mylabel, mytitle)
	
	
	mytitle = "Distribution_of_genes_in_terms_of_fraction_of_retained_introns_"+ cell_type
	rc_cutoff = 0
	get_fraction_retained_intron_per_gene_histogram(reads_on_shared_introns, rc_cutoff, mytitle)
	
	my_intron_dic = get_distribution_ito_tron_number(normalized_rpkms_on_shared_introns)
	my_intron_log_dic = get_distribution_ito_tron_number(log_dist)
	name = "Internal_vs_external_introns_boxplot_" + cell_type
	internal_vs_external_introns_read_densities = Internal_vs_external_trons_boxplot(my_intron_log_dic, name, 4)
	#for number_of_trons in xrange(4,10):
		#generate_characteristic_figures(my_intron_dic, number_of_trons, "intron_read_density_" + str(number_of_trons) + "_" + opt.name)

	return (read_density_relative_fluctuation, internal_vs_external_introns_read_densities)
	

	
			
	
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
	
	# entrez_genes is a EntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
		
	# test module
	if test == 1:
		test_id = 54
		Entrez.test_gene_structure(entrez_gene_collection, test_id)

	
	totalcount_F = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnForwardStrand);
	totalcount_R = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnReverseStrand);
	totalcount = totalcount_F + totalcount_R

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