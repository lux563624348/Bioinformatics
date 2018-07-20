#!/usr/bin/env python

"""
# Calculate read counts on Promoter, Genebody, PromoterGenebody, exonic region at the Transcript level
# For calculation at the gene level, use those under tools

Output: 
name	rc rpkm

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time

import UCSC_revised
import GenomeData
import SeparateByChrom
import associate_tags_with_regions
import Utility_extended
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

def get_read_count_on_genic_regions(geneList, bedFile, fragment_size):
	"""
	only deals with one chrom
	geneList is a UCSC_lite object: name, chrom, strand, txStart, txEnd
	
	Returns three lists: gene name, length,  read count
	"""
	
	(gene_name_list, region_start_list, region_end_list) = get_feature_lists(geneList)	
	tag_position_list = []
	f = open(bedFile,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
	f.close();
	if not Utility_extended.is_list_sorted(tag_position_list):
		tag_position_list.sort()
	#A list, with total tag number on this region, order as the region lists
	read_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, region_start_list, region_end_list)
	
	assert len(gene_name_list) == len(read_count_list)
	region_length_list = [0] * len(gene_name_list)
	for i in xrange(len(gene_name_list)):
		region_length_list[i] = region_end_list[i] - region_start_list[i] 

	return gene_name_list, region_length_list, read_count_list

def test_get_read_count_on_genic_regions(gene_name, gene_name_list, region_length_list, read_count_list):
	list_length = len(gene_name_list)
	for i in xrange(list_length):
		if gene_name_list[i] == gene_name:
			print gene_name_list[i]
			print region_length_list[i]
			print read_count_list[i]
			break
	
	

def get_read_count_on_genes(rawreadfile, fragment_size, knowngenefile, regiontype, promoter_upstream_extension, promoter_downstream_extension):
	"""
	This one provides an integrated module, where the chrom separation etc is done inside. 
	
	Promoter and GeneBody are mutually exclusive. 
	Promoter: TSS-upstreamextention, TSS+downstreamextension
	GeneBody: TSS+downstreamextension, TES
	ExonicRegion: exons of a gene taken together
	PromoterGenebody: Promoter + gene body.
	
	Return: a dictionary with key of gene name and value of read count
	"""
	knowngenes = UCSC_revised.KnownGenes(knowngenefile);
	chroms = knowngenes.keys();
	allowed_region_type = ['Promoter', 'GeneBody', 'PromoterGenebody', 'ExonicRegion'];
	if regiontype not in allowed_region_type:
		print " The allowed region types are Promoter, GeneBody,  PromoterGenebody and ExonicRegion. The region type is not recognized, exiting";
		sys.exit(1);
	if regiontype == 'Promoter':
		region_dic = knowngenes.getPromoters(promoter_upstream_extension, promoter_downstream_extension);
	if regiontype == 'GeneBody':
		region_dic = knowngenes.getGenebodys(promoter_downstream_extension);
	if regiontype == 'PromoterGenebody':
		region_dic = knowngenes.getPromotergenebodys(promoter_upstream_extension);
		
	libName = (rawreadfile).split('/')[-1]
	libName = libName.split('.')[0]
	extension = "-" + libName + ".bed1"
	if Utility_extended.fileExists(rawreadfile):
		SeparateByChrom.separateByChrom(chroms, rawreadfile, extension);
	else:
		print rawreadfile, " not found";
		sys.exit(1)
	
	genes = {}; #for output
	
	for chrom in chroms:
		chrombed = chrom + extension;
		if Utility_extended.fileExists(chrombed):
			gene_coords = knowngenes[chrom]
			if len(gene_coords) > 0:
				if regiontype == 'ExonicRegion':
					(gene_name_list, region_length_list, read_count_list) = get_read_count_on_exons(gene_coords, chrombed, fragment_size)
				else:
					(gene_name_list, region_length_list, read_count_list) = get_read_count_on_genic_regions(region_dic[chrom], chrombed, fragment_size)
				assert len(gene_name_list) == len (region_length_list)
				assert len(gene_name_list) == len (read_count_list)
				#RPKM = [0] * len(gene_name_list)
				for i in xrange(len(gene_name_list)):
					#if region_length_list[i] > 0:
					#	RPKM[i] = read_count_list[i] / (region_length_list[i]/1000.0) / (totalcount/1000000.0)
					genes[gene_name_list[i]] = read_count_list[i]
					#outline = gene_name_list[i] + '\t' + str(read_count_list[i]) + '\t' + str(RPKM[i]) + '\n'
					#f.write(outline)
	
	SeparateByChrom.cleanup(chroms, extension);

	return genes



def get_read_count_on_exons(gene_coords, bedFile, fragment_size):
	"""
	only deals with one chrom
	gene_coords is a list of UCSC object

	Return: three lists: geneName, exonsTotalLength, exonsTotalReadCount
	"""
	tag_position_list = []
	f = open(bedFile,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
	f.close()	
	if not Utility_extended.is_list_sorted(tag_position_list):
		tag_position_list.sort()


	geneName=[]
	exonsTotalLength=[] #used for calculating the RPKM value
	exonsTotalReadCount=[]
	for g in gene_coords:
		geneName.append(g.name)
		if g.exonCount > 0:
			exon_Starts_str = (g.exonStarts.split(','))[:-1] #remove the last '' because the format is '1,2,3,'
			exon_Ends_str = (g.exonEnds.split(','))[:-1]#remove the last '' because the format is '1,2,3,'
			exon_Starts = [int(x) for x in exon_Starts_str]
			exon_Ends = [int(x) for x in exon_Ends_str]
			assert len(exon_Starts) == len(exon_Ends)
			
			totalLength = 0;
			for i in xrange(len(exon_Starts)):
				totalLength +=  exon_Ends[i] - exon_Starts[i]
			exonsTotalLength.append(totalLength)
			
			exon_read_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, exon_Starts, exon_Ends)
			exonsTotalReadCount.append(sum(exon_read_count_list))
		else:
			exonsTotalLength.append(0)
			exonsTotalReadCount.append(0)
	return geneName, exonsTotalLength, exonsTotalReadCount
	
	
	
def main(argv):
	parser = OptionParser()
	
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedFile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'PromoterGenebody' or 'ExonicRegion'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-u", "--promoter_upstream_extension", action="store", type="int", dest="promoter_upstream_extension", help="upstream extension of promoter region from TSS", metavar="<int>")
	parser.add_option("-d", "--promoter_downstream_extension", action="store", type="int", dest="promoter_downstream_extension", help="downstream extension of promoter region from TSS",  metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()

	known_genes = UCSC_revised.KnownGenes(opt.known_genes);
	chroms = known_genes.keys();
	#Promoter and GeneBody are mutually exclusive. 
	#Promoter: TSS-upstreamextention, TSS+downstreamextension
	#GeneBody: TSS+downstreamextension, TES
	#PromoterGenebody: TSS-upstreamextention,  TES.
	allowed_region_type = ['Promoter', 'GeneBody', 'PromoterGenebody', 'ExonicRegion'];
	if opt.region_type not in allowed_region_type:
		print " The allowed region types are Promoter, GeneBody, PromoterGenebody and ExonicRegion. The region type is not recognized, exiting";
		sys.exit(1);
	
	if opt.region_type == 'Promoter':
		region_dic = known_genes.getPromoters(opt.promoter_upstream_extension, opt.promoter_downstream_extension);
	elif opt.region_type == 'GeneBody':
		region_dic = known_genes.getGenebodys(opt.promoter_downstream_extension);
	elif opt.region_type == 'PromoterGenebody':
		region_dic = known_genes.getPromotergenebodys(opt.promoter_upstream_extension);
		
	
	
	libName = (opt.bedFile).split('/')[-1]
	libName = libName.split('.')[0]
	extension = "-" + libName +'.bed1'
	if Utility_extended.fileExists(opt.bedFile):
		SeparateByChrom.separateByChrom(chroms, opt.bedFile,  extension)
	else:
		print opt.bedFile, " not found";
		sys.exit(1)
	
	totalcount = get_total_tag_counts.get_total_tag_counts(opt.bedFile);
	
	f = open(opt.out_file, 'w')
	outline = "# GeneName" +  '\t' + "Read Count" + '\t' + "RPKM" + '\n'
	f.write(outline)
	
	for chrom in chroms:
		chrombed = chrom + extension;
		if Utility_extended.fileExists(chrombed):
			gene_coords = known_genes[chrom]
			if len(gene_coords) > 0:
				if opt.region_type == 'ExonicRegion':
					(gene_name_list, region_length_list, read_count_list) = get_read_count_on_exons(gene_coords, chrombed, opt.fragment_size)
				else:
					(gene_name_list, region_length_list, read_count_list) = get_read_count_on_genic_regions(region_dic[chrom], chrombed, opt.fragment_size)
					#test_get_read_count_on_genic_regions("AAAS", gene_name_list, region_length_list, read_count_list)
					#test_get_read_count_on_genic_regions("AACS", gene_name_list, region_length_list, read_count_list)
				assert len(gene_name_list) == len (region_length_list)
				assert len(gene_name_list) == len (read_count_list)
				RPKM = [0] * len(gene_name_list)
				for i in xrange(len(gene_name_list)):
					if region_length_list[i] > 0:
						RPKM[i] = read_count_list[i] / (region_length_list[i]/1000.0) / (totalcount/1000000.0)
						outline = gene_name_list[i] + '\t' + str(read_count_list[i]) + '\t' + str(RPKM[i]) + '\n'
						f.write(outline)
	f.close()
	
	SeparateByChrom.cleanup(chroms, extension)
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)