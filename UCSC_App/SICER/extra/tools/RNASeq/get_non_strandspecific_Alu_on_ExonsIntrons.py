#!/usr/bin/env python

# Calculate read counts on Promoter, Genebody, PromoterGenebody, exonic region

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import copy
try:
   import cPickle as pickle
except:
   import pickle
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')
sys.path.append('/home/data/hg19/Annotation')

import UCSC_revised
import GenomeData
import SeparateByChrom
import Utility_extended
import RepElements
import Entrez
from Entrez import *

plus = re.compile("\+");
minus = re.compile("\-");


def assign_AluElements_to_intronexons_by_chrom(my_entrez_genes, Alufile_by_chrom, chrom):
	"""
	entrez genes are made sure to be on one chrom, and the bed file are reads for that strand
	The raw read file needs to conform to bed format
	"""
	# Separate by chrom reads
	
	if Utility_extended.fileExists(Alufile_by_chrom) and (chrom in my_entrez_genes.chroms):
		print chrom
		# set up a KnownEntrezGenes Instance for entrez_genes on this particular chrom
		entrez_genes_by_chrom =  Entrez.KnownEntrezGenes([chrom], my_entrez_genes.subset_by_chrom(chrom))
		
		# load in the Alus on chrom
		Alus = RepElements.KnownRepElements.initiate_from_file([chrom], Alufile_by_chrom)
		print "There are %d of elements on %s" %(Alus.number, chrom)
		# Use the mid point of an Alu Element to represent its position, each element is a tuple of (position, id)
		Alu_positions = []
		for myid in Alus.rep_elements.keys():
			element = Alus.rep_elements[myid]
			position = int((element.genoStart + element.genoEnd)/2.0)
			Alu_positions.append((position, myid))
		Alu_positions.sort(key=itemgetter(0))
		
		shared_exon_Alus = {}	#{entrezID:[(region, [Alu_positions])]}
		shared_intron_Alus = {} #{entrezID:[(region, [Alu_positions])]}
		merged_transcript_Alus = {} #{entrezID:[(region, [Alu_positions])]}
		if entrez_genes_by_chrom.num_genes > 0:
			for entrez_id in entrez_genes_by_chrom.entrez_ids:
				
				gene = (entrez_genes_by_chrom.entrez_genes)[entrez_id]
				
				shared_exons = gene.shared_exonic_regions #sorted in absolute coordinate
				if shared_exons == []:# No shared extrons
					shared_exon_Alus[entrez_id] = []
				else:
					# element in result_list has structure: [((start, end), [Alu_positions in range])]
					result_list = Utility_extended.associate_tags_with_regions (Alu_positions, shared_exons) #returns a list, [region, [Alu_positions]]
					shared_exon_Alus[entrez_id] = result_list
					
				shared_introns = gene.shared_intronic_regions #sorted
				if shared_introns == []: # No shared introns
					shared_intron_Alus[entrez_id] = []
				else:
					# element in result_list has structure: [((start, end), [Alus in range])]
					result_list = Utility_extended.associate_tags_with_regions (Alu_positions, shared_introns) #returns a list, [region, [Alu_elements]]
					shared_intron_Alus[entrez_id] = result_list
				
				merged_transcript = gene.boundaries #[(start, end)]
				# element in result_list has structure: [((start, end), [Alus in range])]
				result_list = Utility_extended.associate_tags_with_regions (Alu_positions, merged_transcript)
				merged_transcript_Alus[entrez_id] = result_list
	
	return (shared_intron_Alus, shared_exon_Alus, merged_transcript_Alus)


def test(my_entrez_genes, shared_Alus, id):
	print "Entrez id is %d" % id
	
	gene = my_entrez_genes.entrez_genes[id]
	gene.info()
	
	print
	composite = shared_Alus[id]
	print "The number of shared intronic regions is %d" % len(composite)
	for element in composite:
		region = element[0]		
		Alu_positions = element[1]
		print region, Alu_positions
	


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--AluElementsFile", action="store", type="string", dest="Alus", help="input Alu annotation file for non-strand specific analysis", metavar="<file>")
	parser.add_option("-u", "--entrez_genes_file", action="store", type="string", dest="entrez_collection", metavar="<file>", help="file with curated known genes clustered by entrez ID in pickle format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name for genes and tag numbers")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	startTime = time.time()

	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	# entrez_collection is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_collection, 'rb')
	temp = pickle.load(annotation)
	my_entrez_genes = Entrez.KnownEntrezGenes(chroms, temp) 
	annotation.close()
	
	#test entrez, checks out
	#id = my_entrez_genes.entrez_genes.keys()[0]
	#print id
	#for i in my_entrez_genes.entrez_genes[id].transcripts:
		#print i.getAll()
	
	
	
	lib_name = (opt.Alus).split('/')[-1] # remove directory
	suffix = lib_name.split('.')[-1] # txt
	lib_name = lib_name.split('.')[0] 
	extension = "-" + lib_name +'.' + suffix +"1"
	if Utility_extended.fileExists(opt.Alus):
		if Utility_extended.chrom_files_exist(chroms, extension) != 1:
			# Separate by chrom and sort by start
			print chroms, extension, " files do not exist, separate by chroms. "
			SeparateByChrom.separateByChrom(chroms, opt.Alus, extension)
	else:
		print opt.Alus, " is not found";
		sys.exit(1)

	Alus_in_shared_intron = {}
	Alus_in_shared_exon = {}
	Alus_in_merged_transcript = {}
	
	for chrom in chroms:
		(shared_intron_Alus, shared_exon_Alus, merged_transcript_Alus) = assign_AluElements_to_intronexons_by_chrom(my_entrez_genes, chrom+extension, chrom)
		if chrom == chroms[0]:
			myid = shared_intron_Alus.keys()[0]
			test(my_entrez_genes, shared_intron_Alus, myid)
		Alus_in_shared_intron.update(shared_intron_Alus)
		Alus_in_shared_exon.update(shared_exon_Alus)
		Alus_in_merged_transcript.update(merged_transcript_Alus)
	
	#{entrezID:[(region=(start, end), Alu_count)]}
	Alus_in_shared_intron_dist = {}
	for myid in Alus_in_shared_intron.keys():
		shared_intronic_regions_on_this_gene = Alus_in_shared_intron[myid]
		Alus_on_shared_intronic_regions_on_this_gene = []
		for region in shared_intronic_regions_on_this_gene:
			region_coord, Alu_positions = region
			number_of_Alus = len(Alu_positions)
			Alus_on_shared_intronic_regions_on_this_gene.append((region_coord, number_of_Alus))
		Alus_in_shared_intron_dist[myid] = Alus_on_shared_intronic_regions_on_this_gene
	outname = opt.outfile +  "_Alu_distribution_in_shared_intron.pkl"
	output = open(outname, 'wb')
	pickle.dump(Alus_in_shared_intron_dist, output)
	print "The number of genes output to %s is %d " % (outname, len(Alus_in_shared_intron.keys()))
	output.close()	
	
	#total_intronic_regions = 0
	#for myid in Alus_in_shared_intron.keys():
	#	total_intronic_regions += len(Alus_in_shared_intron[myid])
	#print "There are %d genes with %d shared intronic regions " % (len(Alus_in_shared_intron.keys()),  total_intronic_regions)
	
	#{entrezID:[(region, Alu_positions)]}
	outname = opt.outfile +  "_Alus_in_shared_intron.pkl"
	output = open(outname, 'wb')
	pickle.dump(Alus_in_shared_intron, output)
	print "The number of genes output to %s is %d " % (outname, len(Alus_in_shared_intron.keys()))
	output.close()	
	
	#{entrezID:[(region, Alu_positions)]}
	outname = opt.outfile + "_Alus_in_shared_exon.pkl"
	output = open(outname, 'wb')
	pickle.dump(Alus_in_shared_exon, output)
	print "The number of genes output to %s is %d " % (outname, len(Alus_in_shared_exon.keys()))
	output.close()	
	
	
	
	#Though in this case the structure can be simpler: {entrezID:(region, Alu_count)}, it is better to make the interface uniform.{entrezID:[(region, Alu_count)]}
	Alus_in_merged_transcript_dist = {}
	for myid in Alus_in_merged_transcript.keys():
		assert len(Alus_in_merged_transcript[myid]) == 1
		region = (Alus_in_merged_transcript[myid])[0]
		region_coord, Alu_positions = region
		number_of_Alus = len(Alu_positions)
		Alus_in_merged_transcript_dist[myid] = [(region_coord, number_of_Alus)]
	outname = opt.outfile +  "_Alu_distribution_in_merged_transcript.pkl"
	output = open(outname, 'wb')
	pickle.dump(Alus_in_merged_transcript_dist, output)
	print "The number of genes output to %s is %d " % (outname, len(Alus_in_merged_transcript.keys()))
	output.close()
	
	#{entrezID:[(region, Alu_positions)]}
	outname = opt.outfile + "_Alus_in_merged_transcript.pkl"
	output = open(outname, 'wb')
	pickle.dump(Alus_in_merged_transcript, output)
	print "The number of genes output to %s is %d " % (outname, len(Alus_in_merged_transcript.keys()))
	output.close()	
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)