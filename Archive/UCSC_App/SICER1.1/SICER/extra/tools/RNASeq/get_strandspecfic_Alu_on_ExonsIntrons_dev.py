#!/usr/bin/env python

# Assign Alus on intronic, exonic regions (eg, strand-specific RNA Seq)

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
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/Repeats')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')

import UCSC_revised
import GenomeData
import SeparateByChrom
import Utility_extended
import get_total_tag_counts
import RepElements
import Entrez

plus = re.compile("\+");
minus = re.compile("\-");
		
		

	#SeparateByChrom.cleanup(chroms, rawreadsextension1)
			
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--forwardalufile", action="store", type="string", dest="AlusOnForwardStrand", help="input file for Alus on forward strand", metavar="<file>")
	parser.add_option("-r", "--reversealufile", action="store", type="string", dest="AlusOnReverseStrand", help="input file for Alus on reverse strand", metavar="<file>")
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
	
	# entrez_genes is a EntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_genes = Entrez.EntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
	
	#test module
	Entrez.test(entrez_genes)
	
	
	
	
	
	

	#Clear the file.
	outf = open(opt.outfile, 'w')
	outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	outf.write(outline)
	outf.close()
	
	# The RNA seq data are strand specific. Only use + reads on genes on forward strand, and - reads on genes on reverse strand.
	print "Process genes on forward strand"
	(forward_shared_exon_count, forward_shared_intron_count) = CalculateExonIntrons(entrez_genes_on_forward_strand, opt.ReadsOnForwardStrand, chroms, opt.fragment_size, totalcount, opt.outfile)
	print "Process genes on reverse strand"
	(reverse_shared_exon_count, reverse_shared_intron_count) = CalculateExonIntrons(entrez_genes_on_reverse_strand, opt.AlusOnReverseStrand, chroms, opt.fragment_size, totalcount, opt.outfile)
	
	#combine the densities
	shared_exon_count = {}
	shared_intron_count = {}
	for chrom in chroms:
		# exon
		if chrom in forward_shared_exon_count.keys():
			shared_exon_count[chrom] = forward_shared_exon_count[chrom]
		if chrom in reverse_shared_exon_count.keys():
			shared_exon_count[chrom].update(reverse_shared_exon_count[chrom])	
		# intron
		if chrom in forward_shared_intron_count.keys():
			shared_intron_count[chrom] = forward_shared_intron_count[chrom]
		if chrom in reverse_shared_intron_count.keys():
			shared_intron_count[chrom].update(reverse_shared_intron_count[chrom])
	#store the info in a pickle file
	name = opt.outfile + "_shared_exon_RPKMS.pkl"
	output = open(name, 'wb')
	pickle.dump(shared_exon_count, output)
	output.close()	
	name = opt.outfile + "_shared_intron_RPKMS.pkl"
	output = open(name, 'wb')
	pickle.dump(shared_intron_count, output)
	output.close()
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)