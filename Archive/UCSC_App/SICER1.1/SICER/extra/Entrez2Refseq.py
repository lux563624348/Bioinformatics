#!/usr/bin/env python

"""
Find the refseq ids for a set of given entrez ids

"""
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import bisect
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/Annotation')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')

import Entrez
#import Calculate3UTRUsageIndexFromCuratedGenes
from Entrez import KnownEntrezGenes
from Entrez import EntrezGene 
import gene_set_manipulation

def main(argv):
	parser = OptionParser()
	parser.add_option("-r", "--refseqfile", action="store", type="string", dest="refseq_ucsc_file", help="input ucsc file for annotated genes, eg,  refFlat_hg19_EntrezID.ucsc", metavar="<file>")
	parser.add_option("-i", "--entrezIDfile", action="store", type="string", dest="entrez_ids_file", help="file for entrez ids", metavar="<file>")
	parser.add_option("-u", "--entrez_genes_file", action="store", type="string", dest="entrez_genes", metavar="<file>", help="file with curated known genes clustered by entrez ID in pickle format")
	parser.add_option("-o", "--refseqfile", action="store", type="string", dest="refseq_subset_file", help="ucsc file for refseq transcripts belonging to those entrez_ids", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	# entrez_gene_collection is a KnownEntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()	
	
	entrez_ids = []
	f = open(opt.entrez_ids_file,'r')
	for line in f:
		if not comment.match(line):
			line = line.strip()
			sline = line.split('\t')
			entrez_ids.append(sline[0])
	f.close()
	
	refseq_ids = entrez_gene_collection.get_transcript_ids(entrez_ids)
	gene_set_manipulation.output_UCSCsubset_in_file (opt.refseq_ucsc_file, refseq_ids, opt.refseq_subset_file)


if __name__ == "__main__":
	main(sys.argv)