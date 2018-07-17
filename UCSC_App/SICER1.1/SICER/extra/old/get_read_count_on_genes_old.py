#!/usr/bin/env python

# Calculate read counts on Promoter, Genebody, PromoterGenebody

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


import UCSC
import GenomeData
import SeparateByChrom
import associate_tags_with_regions
import Utility
import get_total_tag_counts

plus = re.compile("\+");
minus = re.compile("\-");

def get_feature_lists(ucsc_list):
	gene_name_list=[];
	region_start_list=[];
	region_end_list=[];
	for item in ucsc_list:
		gene_name_list.append(item.name);
		region_start_list.append(item.txStart);
		region_end_list.append(item.txEnd);
	return gene_name_list, region_start_list, region_end_list

def get_read_count_on_genes(rawreadfile, fragment_size, knowngenefile, regiontype, promoter_upstream_extension, promoter_downstream_extension):
	"""
	Promoter and GeneBody are mutually exclusive. 
	Promoter: TSS-upstreamextention, TSS+downstreamextension
	PromoterGenebody: Promoter + gene body.
	
	Return: a dictionary with key of gene name and value of read count
	"""
	knowngenes = UCSC.KnownGenes(knowngenefile);
	chroms = knowngenes.keys();
	allowed_region_type = ['Promoter', 'GeneBody', 'PromoterGenebody'];
	if regiontype == 'Promoter':
		region_dic = knowngenes.getPromoters(promoter_upstream_extension, promoter_downstream_extension);
	elif regiontype == 'GeneBody':
		region_dic = knowngenes.getGenebodys(promoter_downstream_extension);
	elif regiontype == 'PromoterGenebody':
		region_dic = knowngenes.getPromotergenebodys(promoter_upstream_extension);
	else:
		print " The allowed region types are Promoter, GeneBody and PromoterGenebody. The region type is not recognized, exiting";
		sys.exit(1);
		
	if Utility.fileExists(rawreadfile):
		SeparateByChrom.separateByChrom(chroms, rawreadfile, '.bed1');
	else:
		print rawreadfile, " not found";
		sys.exit(1)
	
	genes = {}; 
	for chrom in chroms:
		(gene_name_list, region_start_list, region_end_list) = get_feature_lists(region_dic[chrom])	
		tag_position_list = []
		read_file = chrom + ".bed1";
		f = open(read_file,'r')
		for line in f:
			if not re.match("#", line):
				line = line.strip()
				sline = line.split()
				tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
		f.close();
		tag_position_list.sort()
		#A list, with total tag number on this region, order as the region lists
		tag_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, region_start_list, region_end_list)
		
		assert len(gene_name_list) == len(tag_count_list)
		for i in range(0, len(gene_name_list)):
			genes[gene_name_list[i]] = tag_count_list[i]
	
	SeparateByChrom.cleanup(chroms,'.bed1');

	return genes
	
	
def main(argv):
	parser = OptionParser()
	
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-f", "--fragmentsize", action="store", type="int", dest="fragment_size", help="fragment size of ChIP-seq reads, in bps", metavar="<int>")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'PromoterGenebody'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-u", "--promoter_upstream_extension", action="store", type="int", dest="promoter_upstream_extension", help="upstream extension of promoter region from TSS", metavar="<int>")
	parser.add_option("-d", "--promoter_downstream_extension", action="store", type="int", dest="promoter_downstream_extension", help="downstream extension of promoter region from TSS",  metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
	
	genes = get_read_count_on_genes(opt.bedfile, opt.fragment_size, opt.known_genes, opt.region_type, opt.promoter_upstream_extension, opt.promoter_downstream_extension)

	totalcount = get_total_tag_counts.get_total_tag_counts(opt.bedfile);

	f = open(opt.out_file, 'w')
	non_zero_genes = 0;
	total_read_count_on_genes = 0; 
	items = genes.items(); # convert to a list
	items.sort(key = operator.itemgetter(1), reverse=True); 
	
	for gene in items:
		if gene[1] > 0: 
			non_zero_genes += 1;
			total_read_count_on_genes += gene[1];
		normalized_count = gene[1]/float(totalcount) * 1000000; 
		f.write(gene[0] + '\t' + str(gene[1]) + '\t' + str(normalized_count) + '\n')
	f.close()
	
	print "Total number of ",  opt.region_type , ": ", len(genes.keys());
	print "Number of ",  opt.region_type , " overlapped with islands: ", non_zero_genes;
	print  total_read_count_on_genes, "of the ", totalcount,  " reads are on ", opt.region_type; 
	
if __name__ == "__main__":
	main(sys.argv)