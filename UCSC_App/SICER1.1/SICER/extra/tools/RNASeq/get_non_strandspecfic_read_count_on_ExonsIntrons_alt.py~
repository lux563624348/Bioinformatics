#!/usr/bin/env python

"""
# Calculate strand_non_specific read counts on introns and exons. 

This is an alternative version. The difference between this version and the other version is that:
In this verson the positive and negative reads are mixed together. The reads after shell sort will not remain sorted after shift because of fragment size. The result is still ok, because there is an additional python sort in the stage of shifting the reads. 

test on 100000 reads gives the same results as the other version. The time is the same too.

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

import GenomeData
import get_total_tag_counts
import Entrez
import get_strandspecific_read_count_on_ExonsIntrons
import Utility_extended

plus = re.compile("\+");
minus = re.compile("\-");
	
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
	
	totalcount = 0
	rawreadslibName1 = (opt.Reads).split('/')[-1]
	rawreadssuffix1 = rawreadslibName1.split('.')[-1] 
	rawreadslibName1 = rawreadslibName1.split('.')[0]
	rawreadsextension1 = "-" + rawreadslibName1 +'.' + rawreadssuffix1 +"1"
	
	if Utility_extended.fileExists(opt.Reads) == 1 :
		totalcount = get_total_tag_counts.get_total_tag_counts(opt.Reads)
	else: # if the all file exist, then use the all file, otherwise use the chrom separated file
		for chrom in chroms:
			chrombed = chrom + rawreadsextension1;
			totalcount1 = get_total_tag_counts.get_total_tag_counts(chrombed);
			print chrom, totalcount1
			totalcount += totalcount1
	

	#Clear the file.
	outf = open(opt.outfile, 'w')
	outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	outf.write(outline)
	outf.close()
	
	(reads_on_shared_exons, reads_on_shared_introns, reads_on_merged_transcripts, summary) = get_strandspecific_read_count_on_ExonsIntrons.CalculateExonIntrons(entrez_gene_collection, opt.Reads, chroms, opt.fragment_size, totalcount, opt.outfile)	
	
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
	
	reads_on_merged_transcripts = {}
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