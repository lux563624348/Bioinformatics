#!/usr/bin/env python

# This module will eventually replace the old get_read_count_on_genes.py
# The reason for the 'GeneBody', 'ExtendedGeneBodys', 'PromoterGenebody' is because of the legacy.


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')


import UCSC_revised
import GenomeData
import SeparateByChrom
import associate_tags_with_regions
import Utility_extended
import get_total_tag_counts
import get_strand_specific_read_count_on_genes

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--readfile", action="store", type="string", dest="ReadFile", help="input bed file for raw reads", metavar="<file>")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--RegionType", action="store", type="string", dest="region_type", metavar="<str>", help=" Region to count tags in: Promoter(txStart-upstream, txStart+downstream), GeneBody (txStart + downstream, txEnd), ExtendedGeneBodys(txStart-upstream, txEnd+downstream), PromoterGenebody(txStart-upstream, txEnd), GeneEnd(txEnd-upstream, txEnd+downstream), ExonicRegion (per exon), IntronicRegion (per intron), Exonictranscript (per transcript), IntronicTranscript (per transcript), 5UTR(txStart, cdsStart), 3UTR(cdsEnd, txEnd)")
	parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension", help="upstream extension of region or location, for Promoter, ExtendedGeneBody, PromoterGenebody and GeneEnd ", metavar="<int>")
	parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension", help="downstream extension of region or location, for Promoter, GeneBody, ExtendedGeneBody and GeneEnd ",  metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
		
	known_genes = UCSC.KnownGenes_list(opt.known_genes);
	chroms = list( set(known_genes.keys()) & set(chroms));
	
	#Promoter and GeneBody are mutually exclusive. 
	#Promoter: TSS-upstreamextention, TSS+downstreamextension
	#GeneBody: TSS-downstreamextension, TES
	#PromoterGenebody: TSS-upstreamextention,TES
	#TES:TES-upstreamextention, TES+downstreamextension
	#
	#ExonicRegions: 
	#IntronicRegions
					
	print "There are ", known_genes.getNumGenes(), " genes. "
	
	totalcount= get_total_tag_counts.get_total_tag_counts(opt.ReadFile);
	
	#Clear the file.
	outf = open(opt.outfile, 'w')
	#outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	#outf.write(outline)
	outf.close()
	
	# The RNA seq data are strand specific. Only use + reads on genes on forward strand, and - reads on genes on reverse strand.
	
	allowed_region_type_1 = ['Promoter', 'GeneBody', 'ExtendedGeneBody', 'PromoterGenebody', 'GeneEnd', 'ExonicRegion', 'IntronicRegion', '5UTR', '3UTR'];
	allowed_region_type_2 = ["ExonicTranscript", "IntronicTranscript"]
	
	if opt.region_type in allowed_region_type_1:
		getReadCount(known_genes, opt.ReadFile, chroms, opt.fragment_size, opt.region_type, opt.upstream_extension, opt.downstream_extension, totalcount, opt.outfile)
	elif opt.region_type in allowed_region_type_2:
		get_read_count_on_onic_transcript(known_genes, opt.ReadFile, chroms, opt.fragment_size, opt.region_type, totalcount, opt.outfile)
	else:
		print " The allowed region types are ", allowed_region_type,  " .The region type is not recognized, exiting"
		sys.exit(1);
	
	print "it took", time.time() - startTime, "seconds."
	
if __name__ == "__main__":
	main(sys.argv)