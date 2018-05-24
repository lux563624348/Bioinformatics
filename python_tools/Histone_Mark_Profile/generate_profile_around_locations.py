#!/usr/bin/env python
import HTSeq
import sys, os
from string import *
from optparse import OptionParser
import numpy

def getTSSProfile(ga, gtf_file, gene_list, window_size, resolution, UpstreamExtension, DownstreamExtension):
	upstream_num_points = UpstreamExtension / resolution
	downstream_num_points = DownstreamExtension / resolution
	total_num_points = upstream_num_points + downstream_num_points + 1
	
	num_transcripts = 0

	tss_pos_set = set()
	profile = numpy.zeros(total_num_points, numpy.int)	
	
	gtffile = HTSeq.GFF_Reader(gtf_file, end_included=True)
	for feature in gtffile:
		if feature.type == "transcript_region" and feature.attr["gene_id"] in gene_list:
			tss_pos_set.add(feature.iv.start_d_as_pos)
			num_transcripts += 1
	
	for tss_pos in tss_pos_set:
		index = 0
		while index < total_num_points:
			count_in_window = 0
			if tss_pos.strand == "+":
				index_pos = tss_pos.pos + (index - upstream_num_points) * resolution
				tss_pos_window_iv = HTSeq.GenomicInterval(tss_pos.chrom, index_pos - window_size / 2, index_pos + window_size / 2)
			elif tss_pos.strand == "-":
				index_pos = tss_pos.pos - (index - upstream_num_points) * resolution
				tss_pos_window_iv = HTSeq.GenomicInterval(tss_pos.chrom, index_pos - window_size / 2 + 1, index_pos + window_size / 2 + 1)
			
			for step_iv, step_count in ga[tss_pos_window_iv].steps():
				count_in_window += step_count * step_iv.length
			profile[index] += count_in_window	
			index += 1
			
	return (num_transcripts, profile)

def getTESProfile(ga, gtf_file, gene_list, window_size, resolution, UpstreamExtension, DownstreamExtension):
	upstream_num_points = UpstreamExtension / resolution
	downstream_num_points = DownstreamExtension / resolution
	total_num_points = upstream_num_points + downstream_num_points + 1
	
	num_transcripts = 0
	
	tes_pos_set = set()
	profile = numpy.zeros(total_num_points, numpy.int)
	
	gtffile = HTSeq.GFF_Reader(gtf_file, end_included=True)
	for feature in gtffile:
		if feature.type == "transcript_region" and feature.attr["gene_id"] in gene_list:
			tes_pos_set.add(feature.iv.end_d_as_pos)
			num_transcripts += 1
	
	for tes_pos in tes_pos_set:
		index = 0
		while index < total_num_points:
			count_in_window = 0
			if tes_pos.strand == "+":
				index_pos = tes_pos.pos + (index - upstream_num_points) * resolution
				tes_pos_window_iv = HTSeq.GenomicInterval(tes_pos.chrom, index_pos - window_size / 2, index_pos + window_size / 2)
			elif tes_pos.strand == "-":
				index_pos = tes_pos.pos - (index - upstream_num_points) * resolution
				tes_pos_window_iv = HTSeq.GenomicInterval(tes_pos.chrom, index_pos - window_size / 2 + 1, index_pos + window_size / 2 + 1)
			
			for step_iv, step_count in ga[tes_pos_window_iv].steps():
				count_in_window += step_count * step_iv.length
			profile[index] += count_in_window	
			index += 1
			
	return (num_transcripts, profile)
	
def main(argv):
	desc="""This is a template for the analysis of aggretated tag distribution with respect to a set of points, such as the TSSs of known genes, with one profile from each strand."""
	parser = OptionParser(description=desc)
	parser.add_option("-g", "--genes_gtf_file", action="store", type="string",
			dest="gtf_file", help="genes gtf file", metavar="<file>")
	parser.add_option("-k", "--known_gene_list", action="store", type="string",
			dest="known_gene_list", help="known gene list", metavar="<file>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bed_file", help="file with tags in bed format", metavar="<file>")
	parser.add_option("-t", "--TypeOfSites", action="store", type="string",
			dest="type", help="TSS, TES, TFBS", metavar="<str>")
	parser.add_option("-o", "--outfile", action="store", type="string",
			dest="outfile", help="outfile name", metavar="<file>")
	parser.add_option("-n", "--normalization", action="store", type="float", dest="norm",
			help="additional normalization in addition to number of sites, number of reads per million and window_size per 1K")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size",
	                  help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-u", "--UpstreamExtension", action="store", type="int",
			dest="upstreamExtension", help="UpstreamExtension", metavar="<int>")
	parser.add_option("-d", "--DownstreamExtension", action="store", type="int",
			dest="downstreamExtension", help="DownstreamExtension", metavar="<int>")
	parser.add_option("-w", "--WindowSize", action="store", type="int",
			dest="window_size", help="window size for averaging. When window size > resolution, there is smoothing", metavar="<int>")	
	parser.add_option("-r", "--resolution", action="store", type="int",
			dest="resolution", help="resolution of the upstream and downstream profile, eg, 5", metavar="<int>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 22:
		parser.print_help()
		sys.exit(1)
		
	gene_list = []	
	f = open(opt.known_gene_list, "r")
	for line in f:
		line = line.strip()
		sline = line.split()
		gene_list.append(sline[0])
		
	fragment_size = opt.fragment_size
	window_size = opt.window_size
	UpstreamExtension = opt.upstreamExtension
	DownstreamExtension = opt.downstreamExtension
	resolution = opt.resolution
	
	print "Upstream extension: %i" % UpstreamExtension
	print "Downstream extension: %i" % DownstreamExtension
	print "Scanning window size: %i" % window_size
	print "Scanning resolution: %i" % resolution
	
	num_tags = 0
	ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
	bedfile = HTSeq.BED_Reader(opt.bed_file)
	for alt in bedfile:
		if alt.iv.strand == "+":
			alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
		elif alt.iv.strand == "-":
			alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
		ga[alt_pos] += 1
		num_tags += 1	
				
	if opt.type == "TSS":	
		(num_transcripts, profile) = getTSSProfile(ga, opt.gtf_file, gene_list, window_size, resolution, UpstreamExtension, DownstreamExtension)
	elif opt.type == "TES":
		(num_transcripts, profile) = getTESProfile(ga, opt.gtf_file, gene_list, window_size, resolution, UpstreamExtension, DownstreamExtension)

			
	normalization = 1.0
	normalization = num_tags/1000000.0
	normalization *= num_transcripts
	normalization *= window_size/1000.0
	normalization *= opt.norm
	
	print "Number of locations: %i" % num_transcripts
	print "Number of reads: %i" % num_tags
	print "Normalization is by total number of reads per million. normalization = %f" % normalization	
	
	f = open(opt.outfile, "w")
	xValues = numpy.arange(-UpstreamExtension / resolution, DownstreamExtension / resolution + 1, resolution)
	normalized_profile = profile / normalization
	for index in range(len(xValues)):
		outline = str(xValues[index]) + "\t" + str(normalized_profile[index]) + "\n"
		f.write(outline)
	f.close()
	
		
if __name__ == "__main__":
	main(sys.argv)

	
	
