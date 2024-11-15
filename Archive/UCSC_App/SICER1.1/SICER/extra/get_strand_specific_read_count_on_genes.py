#!/usr/bin/env python

"""
This module works at the transcript level, not at the entrez gene level
# Calculate strand-specific read counts on Promoter, Genebody, PromoterGenebody, exonic etc region
# The reason for the 'GeneBody', 'ExtendedGeneBodys', 'PromoterGenebody' is because of the legacy.
"""


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

plus = re.compile("\+");
minus = re.compile("\-");



def get_read_count_on_regions(region_list, tag_position_list):
	"""
	only deals with one chrom
	region_list is a list of (start, end)

	Return: TotalLength, TotalReadCount
	"""
	total_length = 0 #used for calculating the RPKM value
	total_readcount = 0
	
	starts = [item[0] for item in region_list]
	ends = [item[1] for item in region_list]
	
	for i in xrange(len(region_list)):
		total_length +=  region_list[i][1] - region_list[i][0] + 1
		
	read_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, starts, ends)
	total_readcount = sum(read_count_list)

	return total_length, total_readcount
	
def get_read_count_list_on_regions(region_list, tag_position_list):
	starts = [item[0] for item in region_list]
	ends = [item[1] for item in region_list]
	lengths = [ (item[1]-item[0]+1) for item in region_list]
	read_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, starts, ends)
	return (lengths, read_count_list)

def getReadCount(KnownGenes, bedfile, chroms, fragment_size, region_type, upstream_extension, downstream_extension, totalcount, out_file):
	"""
	Known genes are made sure to be on one strand, and the bed file are reads for that strand
	The raw read file needs to conform to bed format
	"""
	ReadCount={} # keyed by name, valued by (rc, length, rpkm)
	
	# Separate by chrom reads
	rawreadslibName1 = (bedfile).split('/')[-1]
	rawreadssuffix1 = rawreadslibName1.split('.')[-1] 
	rawreadslibName1 = rawreadslibName1.split('.')[0]
	rawreadsextension1 = "-" + rawreadslibName1 +'.' + rawreadssuffix1 +"1"
	if Utility_extended.fileExists(bedfile):
		if Utility_extended.chrom_files_exist(chroms, rawreadsextension1) != 1:
			# Separate by chrom and sort by start
			print chroms, rawreadsextension1, " files do not exist, separate by chroms and sort each file according to the second column. "
			Utility_extended.separate_by_chrom_sort(chroms, bedfile, rawreadsextension1, [2])
	else:
		print bedfile, " is not found";
		sys.exit(1)

	# dictionary has chrom as key and ucsc_lite object (name, chrom, strand, txStart, txEnd) as values
	if region_type == 'Promoter':
		region_dic = KnownGenes.getPromoters(upstream_extension, downstream_extension)
	elif region_type == 'GeneBody':
		region_dic = KnownGenes.getGenebodys(downstream_extension)
	elif region_type == 'ExtendedGeneBody':
		region_dic = KnownGenes.getExtendedGenebodys(upstream_extension, downstream_extension)
	elif region_type == 'PromoterGenebody':
		region_dic = KnownGenes.getPromotergenebodys(upstream_extension)
	elif region_type == 'GeneEnd':
		region_dic = KnownGenes.getGeneEnds(upstream_extension, downstream_extension)
	elif region_type == 'ExonicRegion':
		region_dic = KnownGenes.getExons()
	elif region_type == 'IntronicRegion':
		region_dic = KnownGenes.getIntrons()
	elif region_type == '5UTR':
		region_dic = KnownGenes.get5UTRs(upstream_extension, downstream_extension)
	elif region_type == '3UTR':
		region_dic = KnownGenes.get3UTRs(upstream_extension, downstream_extension)
	else:
		print region_type, "is not recognized"
		exit(1)
		
	outf = open(out_file, 'a')
	
	for chrom in chroms:
		chrombed = chrom + rawreadsextension1;
		if Utility_extended.fileExists(chrombed) and (chrom in KnownGenes.keys()):
			tag_position_list = []
			inf = open(chrombed,'r')
			for line in inf:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
			inf.close()	
			if Utility_extended.is_list_sorted(tag_position_list) != 1:
				tag_position_list.sort()
			
			if len(region_dic[chrom]) > 0:
				for region in region_dic[chrom]:
					thisregion = [(region.txStart, region.txEnd)]
					(total_length, rc) = get_read_count_on_regions(thisregion, tag_position_list)
					if total_length > 0:
						RPKM = rc * (1000.0/total_length) * (1000000/float(totalcount))
					else:
						assert rc<0.01
						RPKM = 0
					outline = str(region.name)  + '\t' + str(rc) + '\t' + str(total_length) + '\t' + str(RPKM) + '\n'
					outf.write(outline)
					ReadCount[region.name] = (rc, total_length, RPKM)
	outf.close()
		
	#SeparateByChrom.cleanup(chroms, rawreadsextension1)
	
	return ReadCount



def get_read_count_on_onic_transcript(KnownGenes, bedfile, chroms,  fragment_size, region_type,  totalcount, out_file):
	"""
	Return: a dictionary keyed by geneName valued by TotalReadCount,TotalLength, RPKM
	"""
	
	ReadCount={} # keyed by name, valued by (rc, length, rpkm)
	
	
	# Separate by chrom reads
	rawreadslibName1 = (bedfile).split('/')[-1]
	rawreadssuffix1 = rawreadslibName1.split('.')[-1] 
	rawreadslibName1 = rawreadslibName1.split('.')[0]
	rawreadsextension1 = "-" + rawreadslibName1 +'.' + rawreadssuffix1 +"1"
	if Utility_extended.fileExists(bedfile):
		if Utility_extended.chrom_files_exist(chroms, rawreadsextension1) != 1:
			# Separate by chrom and sort by start
			print chroms, rawreadsextension1, " files do not exist, separate by chroms and sort each file according to the second column. "
			Utility_extended.separate_by_chrom_sort(chroms, bedfile, rawreadsextension1, [2]) # sort by start
	else:
		print bedfile, " is not found";
		sys.exit(1)
	
	outf = open(out_file, 'a')
	for chrom in chroms:
		chrombed = chrom + rawreadsextension1;
		if Utility_extended.fileExists(chrombed) and (chrom in KnownGenes.keys()):
			tag_position_list = []
			inf = open(chrombed,'r')
			for line in inf:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_position_list.append(associate_tags_with_regions.tag_position(sline, fragment_size))
			inf.close()	
			if Utility_extended.is_list_sorted(tag_position_list) != 1:
				tag_position_list.sort()
			
			for gene in KnownGenes[chrom]:
				if region_type == "ExonicTranscript":
					ons = gene.getExons()
				elif region_type == "IntronicTranscript":	
					ons = gene.getIntrons()
				else:
					print region_type, "is not recognized."
					exit(1)
				if len(ons>0):
					(total_length, rc) = get_read_count_on_regions(ons, tag_position_list)
					RPKM = rc * (1000.0/total_length) * (1000000/float(totalcount))
				else:
					total_length = 0 
					rc = 0
					RPKM = 0
				outline = str(gene.name)  + '\t' + str(rc) + '\t' + str(total_length) + '\t' + str(RPKM) + '\n'
				outf.write(outline)
				ReadCount[region.name] = (rc, total_length, RPKM)
	outf.close()
	
	#SeparateByChrom.cleanup(chroms, rawreadsextension1)
	
	return ReadCount	
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--forwardreadfile", action="store", type="string", dest="ReadsOnForwardStrand", help="input bed file for RNASeq raw reads on forward strand", metavar="<file>")
	parser.add_option("-n", "--reversereadfile", action="store", type="string", dest="ReadsOnReverseStrand", help="input bed file for RNASeq raw reads on reverse strand", metavar="<file>")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--RegionType", action="store", type="string", dest="region_type", metavar="<str>", help=" Region to count tags in: Promoter(txStart-upstream, txStart+downstream), GeneBody (txStart + downstream, txEnd), ExtendedGeneBodys(txStart-upstream, txEnd+downstream), PromoterGenebody(txStart-upstream, txEnd), GeneEnd(txEnd-upstream, txEnd+downstream), ExonicRegion (per exon), IntronicRegion (per intron), Exonictranscript (per transcript), IntronicTranscript (per transcript), 5UTR(txStart, cdsStart), 3UTR(cdsEnd, txEnd)")
	parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension", help="upstream extension of region or location, for Promoter, ExtededGeneBody, PromoterGenebody and GeneEnd ", metavar="<int>")
	parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension", help="downstream extension of region or location, for Promoter, GeneBody, ExtendedGeneBody and GeneEnd ",  metavar="<int>")
	parser.add_option("-o", "--out_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 18:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
		
	known_genes = UCSC_revised.KnownGenes(opt.known_genes);
	chroms = list( set(known_genes.keys()) & set(chroms));
	
	#Promoter and GeneBody are mutually exclusive. 
	#Promoter: TSS-upstreamextention, TSS+downstreamextension
	#GeneBody: TSS-downstreamextension, TES
	#PromoterGenebody: TSS-upstreamextention,TES
	#TES:TES-upstreamextention, TES+downstreamextension
	#ExonicRegions: 
	#IntronicRegions

	total_forward = 0
	total_reverse = 0

	known_genes.output( known_genes.get_strand_specific_genes("+"), "temp_p_genes.ucsc")
	genes_on_forward_strand = UCSC_revised.KnownGenes("temp_p_genes.ucsc")
	print "There are %d genes on forward strand. " %genes_on_forward_strand.getNumGenes()
	
	known_genes.output( known_genes.get_strand_specific_genes("-"), "temp_n_genes.ucsc")
	genes_on_reverse_strand = UCSC_revised.KnownGenes("temp_n_genes.ucsc")
	print "There are %d genes on reverse strand. " %genes_on_reverse_strand.getNumGenes()
	
	totalcount_F = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnForwardStrand);
	totalcount_R = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnReverseStrand);
	totalcount = totalcount_F + totalcount_R
	
	#Clear the file.
	outf = open(opt.out_file, 'w')
	#outline = "# Entrez ID \t Merged Exon Read Count \t Merged Exon Length \t Merged Exon RPKM \t Shared Exon Read Count \t  Shared Exon Length \t Shared Exon RPKM \t Shared Intron Read Count \t Share Intron Length \t Shared Intron RPKM \t Merged Transcript Read Count \t Merged Transcript Length \t Merged Transcript RPKM \t RefSeq IDs \t Gene Symbols \n"
	#outf.write(outline)
	outf.close()
	
	# The RNA seq data are strand specific. Only use + reads on genes on forward strand, and - reads on genes on reverse strand.
	
	allowed_region_type_1 = ['Promoter', 'GeneBody', 'ExtendedGeneBody', 'PromoterGenebody', 'GeneEnd', 'ExonicRegion', 'IntronicRegion', '5UTR', '3UTR'];
	allowed_region_type_2 = ["ExonicTranscript", "IntronicTranscript"]
	
	if opt.region_type in allowed_region_type_1:
		print "Process genes on forward strand"
		getReadCount(genes_on_forward_strand, opt.ReadsOnForwardStrand, chroms, opt.fragment_size, opt.region_type, opt.upstream_extension, opt.downstream_extension, totalcount, opt.out_file)
		print "Process genes on reverse strand"
		getReadCount(genes_on_reverse_strand, opt.ReadsOnReverseStrand, chroms, opt.fragment_size, opt.region_type, opt.upstream_extension, opt.downstream_extension, totalcount, opt.out_file)
	elif opt.region_type in allowed_region_type_2:
		print "Process genes on forward strand"
		get_read_count_on_onic_transcript(genes_on_forward_strand, opt.ReadsOnForwardStrand, chroms, opt.fragment_size, opt.region_type, totalcount, opt.out_file)
		print "Process genes on reverse strand"
		get_read_count_on_onic_transcript(genes_on_reverse_strand, opt.ReadsOnReverseStrand, chroms, opt.fragment_size, opt.region_type, totalcount, opt.out_file)
	else:
		print " The allowed region types are ", allowed_region_type,  ". The region type is not recognized, exiting"
		sys.exit(1);
	
	print "it took", time.time() - startTime, "seconds."
	
if __name__ == "__main__":
	main(sys.argv)