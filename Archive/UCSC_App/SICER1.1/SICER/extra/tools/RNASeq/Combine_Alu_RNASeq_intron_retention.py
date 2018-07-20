#!/usr/bin/env python

""" 
The starting point for RNA-Seq is obtained from get_read_count_on_ExonsIntrons.py: rpkm values for all shared introns and exons

THe starting point for Alu is obtained from get_non_strandspecific_Alu_on_ExonsIntrons_dev.py
"""
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time


import numpy
import matplotlib.pyplot as plt
import matplotlib


import rpy2.robjects as robjects
r = robjects.r
from rpy2.robjects.packages import importr


try:
   import cPickle as pickle
except:
   import pickle
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/Annotation')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')
sys.path.append('/home/data/hg19/JunZhu/modules')

import GenomeData
import Utility_extended
import AnalyzeRNASeqPASeq
from  Entrez import *
import Entrez

plus = re.compile("\+");
minus = re.compile("\-");

def test(read_density_distribution, my_entrez_id):
	print (read_density_distribution)[my_entrez_id]


def readin_Alu(mypickle):
	#The data structure of the returned object is {entrez_id:[((start,end), alu_count)]}
	distribution = open(mypickle, 'rb')
	alu_density_distribution = pickle.load(distribution) 
	print "There are ", len(alu_density_distribution.keys()), " genes in ", mypickle
	return alu_density_distribution

def find_genic_Alu_density(Alus_distribution_in_merged_transcript):
	"""
	#Alus_distribution_in_merged_transcript: {entrezID:(region, Alu_count)}, but only one region, 
	"""
	genic_Alu_density = {}
	genic_Alu_count = {}
	for myid in Alus_distribution_in_merged_transcript.keys():
		on_this_gene = Alus_distribution_in_merged_transcript[myid]
		region, Alu_count = on_this_gene
		#print type(region), type(Alu_count)
		
		start, end = region
		#print type(start), type(end)
		length = end - start
		
		density = float(Alu_count)/length
		genic_Alu_density[myid] = density
		genic_Alu_count[myid] = Alu_count
	return genic_Alu_count, genic_Alu_density
		
def readin_RNASeq(mypickle):
	#The data structure from the pickle is  {chrom:{entrez_id:[rpkms]}}
	#The data structure of the returned object is {entrez_id:[rpkms]}
	distribution = open(mypickle, 'rb')
	read_density_distribution = pickle.load(distribution) 
	distribution.close()
	flattened = {}
	#print "The keys of " + opt.RNA_Seq_intron_pickle + " are: ", read_density_distribution.keys()
	total = 0
	for chrom in read_density_distribution.keys():
		flattened.update(read_density_distribution[chrom])
		total += len(read_density_distribution[chrom].keys())
	
	assert (total == len(flattened.keys()))
	
	#for myid in flattened.keys()[:10]:
		#print flattened[myid]
		
	
	print "There are ", total, " genes in ", mypickle
	return flattened

def find_genic_iri(intron_read_density_distribution, exon_read_density_distribution, idset, exon_read_density_cutoff=0.15):
	total = 0
	genic_iri = {} #{entrez_id:average_rpkm}
	id_set = list( set((intron_read_density_distribution).keys()) & set((exon_read_density_distribution).keys()) & set(idset))
	for myid in id_set:
		exon_rpkms = (exon_read_density_distribution)[myid]
		intron_rpkms = (intron_read_density_distribution)[myid] 
		average_exon_rpkm = numpy.average(exon_rpkms)
		if average_exon_rpkm > exon_read_density_cutoff:
			total += 1
			average_intron_rpkm = numpy.average(intron_rpkms)
			genic_iri[myid] = average_intron_rpkm/float(average_exon_rpkm)
	print "There are ", total, "entrez_ids whose expression >= ", exon_read_density_cutoff, " RPKM"
	return genic_iri

def find_genic_iri_distribution(intron_read_density_distribution, exon_read_density_distribution, idset,  exon_read_density_cutoff=0.15):
	"""
	normalize the intron read densities by the average exon read density of that gene
	exclude genes whose expression is below cutoff.
	"""
	total = 0
	genic_iris = {} #{entrez_id:[rpkms]}
	id_set = list( set((intron_read_density_distribution).keys()) & set((exon_read_density_distribution).keys()) & set(idset))
	for myid in id_set:
		exon_rpkms = (exon_read_density_distribution)[entrez_id]
		average_exon_rpkm = numpy.average(exon_rpkms)
		if average_exon_rpkm >= exon_read_density_cutoff:
			total += 1
			intron_rpkms = (intron_read_density_distribution)[entrez_id]
			genic_iris[entrez_id] = [rpkm/average_exon_rpkm for rpkm in intron_rpkms]
	print "There are ", total, "entrez_ids whose expression >= ", exon_read_density_cutoff, " RPKM"
	return genic_iris

def log_transform_read_density_distribution(rdd, pc=0.0000000001):
	log_read_density_distribution = {}
	for entrez_id in rdd.keys():
		rpkms = rdd[entrez_id]
		log_read_density_distribution[entrez_id] = [log(item+pc,2) for item in rpkms]
	return log_read_density_distribution
			
def find_simple_gene_IRI_distribution(intron_read_density_distribution, exon_read_density_distribution, simple_gene_list, exon_read_density_cutoff=0.15, pc = 0.000000001):
	"""
	simple_gene_list: Simple genes are genes whose entrez id has only on refseq id
	
	For these genes, the introns are in between exons.
	
	returns a {entrez_id:[iri]}
	"""
	total = 0
	iri = {}
	mykeys = list (set((intron_read_density_distribution).keys()) & set (simple_gene_list) )
	for entrez_id in mykeys:
		intron_rpkms = (intron_read_density_distribution)[entrez_id]
		exon_rpkms = (exon_read_density_distribution)[entrez_id]
		if len(exon_rpkms) != len(intron_rpkms) + 1:
			print entrez_id, len(exon_rpkms), len(intron_rpkms)
			print exon_rpkms
			print intron_rpkms
		assert (len(exon_rpkms) == len(intron_rpkms) + 1)
		if numpy.average(exon_rpkms) >= exon_read_density_cutoff:
			total += 1
			iri[entrez_id] = []
			for i in xrange(len(intron_rpkms)):
				if intron_rpkms[i] == 0.0:
					myiri = 0
				else:
					neighbor_exon_rpkm = (exon_rpkms[i] + exon_rpkms[i+1])/2.0
					myiri = (intron_rpkms[i] + pc) / (neighbor_exon_rpkm + pc)
				iri[entrez_id].append(myiri)
	print "There are ", total, " simple genes whose expression >= ", exon_read_density_cutoff, " RPKM."
	return iri

def combine_iri_Alu_genic_level(intron_read_density_distribution, exon_read_density_distribution, Alus_distribution_in_merged_transcript, id_subset, exon_read_density_cutoff=0.15):

	genic_iri = find_genic_iri(intron_read_density_distribution, exon_read_density_distribution, id_subset,exon_read_density_cutoff)
	genic_Alu_count, genic_Alu_density = find_genic_Alu_density(Alus_distribution_in_merged_transcript)
	
	idset = list( set(genic_iri.keys()) & set(genic_Alu_count.keys()) & set(id_subset))
	
	genic_iri_list = []
	genic_Alu_count_list = []
	genic_Alu_density_list = []
	for myid in idset:
		genic_iri_list.append(genic_iri[myid])
		genic_Alu_count_list.append(genic_Alu_count[myid])
		genic_Alu_density_list.append(genic_Alu_density[myid])
	AnalyzeRNASeqPASeq.scatterplot(genic_iri_list, genic_Alu_count_list, "genic_iri_vs_Alu_count", xscale='log', yscale='log')
	AnalyzeRNASeqPASeq.scatterplot(genic_iri_list, genic_Alu_density_list, "genic_iri_vs_Alu_density", xscale='log', yscale='log')

def combine_iri_Alu_intron_level(intron_read_density_distribution, exon_read_density_distribution, alu_density_distribution, id_subset, exon_read_density_cutoff=0.15, pc = 0.000000001):			
	"""
	Collect all the introns, and draw a scatter plot for alu_count vs iri
	alu_density_distribution:{entrez_id:[((start,end), alu_count)]}
	"""
	#{entrez_id:[iri]}
	idset = list(set(alu_density_distribution.keys())  & set(id_subset))
	my_iris = find_simple_gene_IRI_distribution(intron_read_density_distribution, exon_read_density_distribution, idset,  exon_read_density_cutoff, pc)
	
	idset = list( set(my_iris.keys()) & set(alu_density_distribution.keys())  & set(id_subset))
	all_introns_iri = []
	all_introns_Alus = []
	for myid in idset:
		assert len(my_iris[myid]) == len(alu_density_distribution[myid])
		for intron_index in xrange(len(my_iris[myid])):
			all_introns_iri.append(my_iris[myid][intron_index])
			alu_on_this_intron = alu_density_distribution[myid][intron_index]
			all_introns_Alus.append(alu_on_this_intron[1])
	AnalyzeRNASeqPASeq.scatterplot(all_introns_iri, all_introns_Alus, "intronic_iri_vs_Alu_count", xscale='log', yscale='log')
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-u", "--annotation_pickle_file", action="store", type="string", dest="annotation", metavar="<file>", help="annotation for strand information")
	parser.add_option("-a", "--rnaseq_intron_pickle_file", action="store", type="string", dest="RNA_Seq_intron_pickle", metavar="<file>", help="read densities for individual shared intronic regions in pickle format")
	parser.add_option("-b", "--rnaseq_exon_pickle_file", action="store", type="string", dest="RNA_Seq_exon_pickle", metavar="<file>", help="read densities for individual shared exonic regions in pickle format")
	parser.add_option("-d", "--alu_distribution_intron_pkl", action="store", type="string", dest="alu_distribution_intron_pkl", metavar="<file>", help="Alu densities for individual shared intronic regions in pickle format")
	parser.add_option("-e", "--alu_distribution_in_merged_transcript_pkl", action="store", type="string", dest="alu_distribution_in_merged_transcript_pkl", metavar="<file>", help="Alu counts for transcript regions in pickle format")
	parser.add_option("-f", "--ids", action="store", type="string", dest="id_subset_file", metavar="<file>", help="file that records ids of interest", default="")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
		parser.print_help()
		sys.exit(1)
        		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);   	
        	
	print "Read in RNASeq information"
	# read in from intron pickle resulted in {entrezID:[rpkms]}
	intron_read_density_distribution =  readin_RNASeq(opt.RNA_Seq_intron_pickle)
	# read in from exon pickle resulted in {entrezID:[rpkms]}
	exon_read_density_distribution =  readin_RNASeq(opt.RNA_Seq_exon_pickle)
	
	# read in annotation file for strand and gene cluster information
	print "Read in annotation"
	annotation = open(opt.annotation, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
	
	# read in the Alus
	print "Read in Alu information" 
	Alu_distribution_on_introns =  readin_Alu(opt.alu_distribution_intron_pkl)
	Alu_distribution_on_transcripts = readin_Alu(opt.alu_distribution_in_merged_transcript_pkl)
	
	print "Load in ids"
	if opt.id_subset_file != "":
		id_set = []
		f = open(opt.id_subset_file,'r')
		for line in f:
			if not comment.match(line):
				line = line.strip()
				sline = line.split('\t')
				id_set.append( int(sline[0]))
		f.close()
		print "There are %d ids in %s" % (len(id_set), opt.id_subset_file)
	

	pc = 0.000000001
	exon_read_density_cutoff = 1
	combine_iri_Alu_genic_level(intron_read_density_distribution, exon_read_density_distribution, Alu_distribution_on_transcripts, id_set, exon_read_density_cutoff)
	
	combine_iri_Alu_intron_level(intron_read_density_distribution, exon_read_density_distribution, Alu_distribution_on_introns, id_set, exon_read_density_cutoff, pc)
	
	
	
if __name__ == "__main__":
	main(sys.argv)