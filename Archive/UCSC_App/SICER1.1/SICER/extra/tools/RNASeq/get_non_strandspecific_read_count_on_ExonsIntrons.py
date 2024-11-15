#!/usr/bin/env python
"""

Calculate strand non-specific read counts on intronic, exonic regions with strand non-specific reads (eg, strand-non-specific RNA Seq). This code is based on the strand-specific version.

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
try:
   import cPickle as pickle
except:
   import pickle
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/Annotation')

import UCSC_revised
import GenomeData
import SeparateByChrom
import Utility_extended
import get_total_tag_counts
import get_strandspecific_read_count_on_ExonsIntrons
import Entrez
from Entrez import EntrezGene
from Entrez import KnownEntrezGenes

plus = re.compile("\+");
minus = re.compile("\-");


def calculate_non_strandspecific_rc_on_ExonIntrons(entrez_genes, bedfile, chroms,  fragment_size):
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
		p_file_name = bedfile + "_P"
		n_file_name = bedfile + "_n"
		Utility_extended.separate_by_strand(bedfile, p_file_name, n_file_name) #partition the bed file into reads in positive strand and negative strand
		##################################################################3
		#The column numbers are 1 based instead of 0 based!
		#For positive strand
		start_index_P = 2
		#For negative strand
		start_index_N = 3
		##################################################################3
		p_totalcount = get_total_tag_counts.get_total_tag_counts(p_file_name)
		(forward_reads_on_shared_exons, forward_reads_on_shared_introns, forward_reads_on_merged_transcripts, forward_summary) = get_strandspecific_read_count_on_ExonsIntrons.calculateExonIntrons(entrez_genes, p_file_name, start_index_P, chroms,  fragment_size, p_totalcount, None)
		
		n_totalcount = get_total_tag_counts.get_total_tag_counts(n_file_name)
		(reverse_reads_on_shared_exons, reverse_reads_on_shared_introns, reverse_reads_on_merged_transcripts, reverse_summary) = get_strandspecific_read_count_on_ExonsIntrons.calculateExonIntrons(entrez_genes, n_file_name, start_index_N, chroms,  fragment_size, n_totalcount, None)
		
		all_reads_on_shared_exons = {} # {entrezID:[((start, end), read_count)]}
		all_reads_on_shared_introns = {} # {entrezID:[((start, end), read_count)]}
		all_reads_on_merged_transcripts = {} #{entrezID:[((start, end), read_count)]}
		all_summary = {} # {entrezID:{attributes}}
	
		all_reads_on_shared_exons = combine_rc(forward_reads_on_shared_exons, reverse_reads_on_shared_exons)
		all_reads_on_shared_introns = combine_rc(forward_reads_on_shared_introns, reverse_reads_on_shared_introns)
		all_reads_on_merged_transcripts = combine_rc(forward_reads_on_merged_transcripts, reverse_reads_on_merged_transcripts)
		all_summary = combine_summary(forward_summary, reverse_summary, p_totalcount, n_totalcount)
		
	SeparateByChrom.cleanup(chroms, extension)
	return (all_reads_on_shared_exons, all_reads_on_shared_introns, all_reads_on_merged_transcripts, all_summary)

def combine_rc(forward, reverse):
	"""
	forward:  {entrezID:[((start, end), read_count)]}
	reverse: {entrezID:[((start, end), read_count)]}
	"""
	print "There are %d entrez ids in forward" %len(forward.keys()) 
	print "There are %d entrez ids in reverse" %len(reverse.keys()) 
	myentrez_ids = list(set(forward.keys()) & set(reverse.keys()))
	print "The overlap between forward and reverse has %d entrez ids" %len(myentrez_ids)
		
	all_reads = {}
	for myid in myentrez_ids:
		all_reads[myid] = []
		for i in xrange(len(forward[myid])):
			myrange = (forward[myid])[i][0]
			rc = forward[myid][i][1]  + reverse[myid][i][1]
			all_reads[myid].append((myrange, rc))
	return all_reads
	
def combine_summary(forward, reverse, p_totalcount, n_totalcount):
	"""
	forward:  {entrezID:{attribute:value}}
	reverse: {entrezID:{attribute:value}}
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
	print "There are %d entrez ids in forward_summary" %len(forward.keys()) 
	print "There are %d entrez ids in reverse_summary" %len(reverse.keys()) 
	myentrez_ids = list(set(forward.keys()) & set(reverse.keys()))
	print "The overlap between forward_summary and reverse summary has %d entrez ids" %len(myentrez_ids)
		
	all_summary = {}
	for myid in myentrez_ids:
		all_summary[myid] = {}
		for feature in forward[myid].keys():
			all_summary[myid][feature] = forward[myid][feature]
		for feature in ["merged_exons_rc", "shared_exons_rc", "shared_introns_rc",  "merged_transcript_rc"]:
			all_summary[myid][feature] = forward[myid][feature] + reverse[myid][feature]
		for feature in ["merged_exon_RPKM", "shared_exon_RPKM", "shared_intron_RPKM", "merged_transcript_RPKM"]:
			all_summary[myid][feature] = (forward[myid][feature] * p_totalcount + reverse[myid][feature] * n_totalcount )/(p_totalcount*1.0 + n_totalcount)
	return all_summary
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-r", "--readfile", action="store", type="string", dest="Reads", help="input bed file for non-strand specific raw reads", metavar="<file>")
	parser.add_option("-g", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determines the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-u", "--entrez_genes_file", action="store", type="string", dest="entrez_genes", metavar="<file>", help="file with curated known genes clustered by entrez ID in pickle format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name for genes and tag numbers")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	startTime = time.time()

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
		test_id = 54
		Entrez.test_gene_structure(entrez_gene_collection, test_id)
	
	
	rawreadslibName1 = (opt.Reads).split('/')[-1]
	rawreadssuffix1 = rawreadslibName1.split('.')[-1] 
	rawreadslibName1 = rawreadslibName1.split('.')[0]
	rawreadsextension1 = "-" + rawreadslibName1 +'.' + rawreadssuffix1 +"1"
	
	totalcount = 0
	if Utility_extended.fileExists(opt.Reads) == 1 :
		totalcount = get_total_tag_counts.get_total_tag_counts(opt.Reads)
	else: # if the all file exist, then use the all file, otherwise use the chrom separated file
		for chrom in chroms:
			chrombed = chrom + rawreadsextension1;
			totalcount1 = get_total_tag_counts.get_total_tag_counts(chrombed);
			print chrom, totalcount1
			totalcount += totalcount1
	
	
	(reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary) = calculate_non_strandspecific_rc_on_ExonIntrons(entrez_gene_collection, opt.Reads, chroms,  opt.fragment_size) 
	
	
	#Clear the file.
	outf = open(opt.outfile, 'w')
	outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	outf.write(outline)
	for entrez_id in entrez_gene_collection.entrez_ids:
		gene = (entrez_gene_collection.entrez_genes)[entrez_id]
		gene_symbol = []
		for transcript in gene.transcripts:
			if transcript.additional_annotations[0] not in gene_symbol:
				gene_symbol.append(transcript.additional_annotations[0])
		outline = str(entrez_id) + '\t' + str(summary[entrez_id]["merged_exons_rc"]) + '\t' + str(summary[entrez_id]["merged_exons_total_length"]) + '\t' + str(summary[entrez_id]["merged_exon_RPKM"]) + '\t' + str(summary[entrez_id]["shared_exons_rc"]) + '\t' + str(summary[entrez_id]["shared_exons_total_length"]) + '\t' + str(summary[entrez_id]["shared_exon_RPKM"]) + '\t' + str(summary[entrez_id]["shared_introns_rc"]) + '\t' + str(summary[entrez_id]["shared_introns_total_length"]) + '\t' +str(summary[entrez_id]["shared_intron_RPKM"]) + '\t' + str(summary[entrez_id]["merged_transcript_rc"]) + '\t' + str(summary[entrez_id]["merged_transcript_length"]) + '\t' + str(summary[entrez_id]["merged_transcript_RPKM"]) + '\t' + ','.join([transcript.name for transcript in gene.transcripts]) + '\t' + ','.join(gene_symbol) + '\n'
		outf.write(outline)
	outf.close()
	
	# {entrezID:[((start, end), read_count)]}
	name = opt.outfile + "_shared_exons.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_shared_exons, output)
	output.close()
	
	# {entrezID:[((start, end), read_count)]}
	name = opt.outfile + "_shared_introns.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_shared_introns, output)
	output.close()
	
	#store the info in a pickle file
	name = opt.outfile + "_merged_transcripts.pkl"
	output = open(name, 'wb')
	pickle.dump(reads_on_merged_transcripts, output)
	output.close()
	
	name = opt.outfile + "_summary.pkl"
	output = open(name, 'wb')
	pickle.dump(summary, output)
	output.close()
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)
