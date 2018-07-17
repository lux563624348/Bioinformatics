#!/usr/bin/env python

"""
03/2011

This is a template for the analysis of tag distribution with respect
to a set of points, such as the TSSs of known genes.

desired functionality: combination of 5' profile and 3' profile 
"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import bisect
import time

import Utility
import BED;
import UCSC;
import GenomeData 
import SeparateByChrom
import associate_tags_with_regions


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

# The strand orientation is taken care of in BED2, which is assumed for bed_vals.
def breakUpStrands(bed_vals):
    plus_starts = [];
    minus_starts = [];
    for b in bed_vals:
        if plus.match(b.strand):
            plus_starts.append(b.start);
        elif minus.match(b.strand):
            minus_starts.append(b.start);
    #print len(plus_starts), len(minus_starts)
    return (plus_starts, minus_starts)


def getProfileNearPosition(position, orientation, upstream_length, downstream_length, resolution, window_size, tag_positions):
	"""
	Draw profiles around a set of positions defined on a particular chrom
	
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	position is a position with orientation
	tag_positions is a list of positions on the same chrom. 
		Tags also have orientation, they are taken into account elsewhere.
	Output is not normalized at all
	
	"""
	length = upstream_length + downstream_length;
	numPoints = length/resolution;
	half_partition = int(resolution/2.0)
	scores = []
	for i in range(numPoints):
		scores.append(0.0);
	if len(tag_positions) > 0:
		if plus.match(orientation): 
			for index in range(numPoints):
				# xValue is the mid point of each partition
				xValue = position - upstream_length + half_partition + index * resolution;
				start = xValue - int(window_size/2.0);
				end = start + window_size - 1;
				count = associate_tags_with_regions.countTagsInWindow(start, end, tag_positions)
				scores[index] += count
				#score[index] += float(count)/float(window_size)
		elif minus.match(orientation): 
			for index in range(numPoints):
				# xValue is the mid point of each partition
				xValue = position + upstream_length - half_partition - index * resolution;
				end = xValue + int(window_size/2.0);
				start = end - window_size + 1;
				count = associate_tags_with_regions.countTagsInWindow(start, end, tag_positions)
				scores[index] += count
		else:
			print "Wrong value for orientation, which can only be + or -";
			sys.exit (1);
	return scores;
	
	
	

def getTSSProfile(coords, upstream_length, downstream_length, resolution, window_size, bed_vals):
	"""
	Draw profiles around a set of TSSes 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary of UCSC-type objects keyed by chrom
	bed-vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is pure count, not normalized at all
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		exit (1)
	numPoints = length/resolution;
	plus_scores = [];
	minus_scores = [];
    	for i in range(int(numPoints)):
        	plus_scores.append(0);
        	minus_scores.append(0);
	
	for chrom in coords.keys():
        	genes = coords[chrom];
        	if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = breakUpStrands(bed_vals[chrom]);
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			for g in genes:
				if plus.match(g.strand):
					temp_score = getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)
					for index in xrange(int(numPoints)):
						plus_scores[index] += temp_score[index];
					temp_score = getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)	
					for index in xrange(int(numPoints)):
						minus_scores[index] += temp_score[index];	
				elif minus.match(g.strand):
					## notice, swapped strands for tags because we're working on the crick strand now
					# ie, for genes on the minus strand, plus reads --> minus reads, minus-reads --> plus reads
					temp_score = getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)
					for index in xrange(int(numPoints)):
						plus_scores[index] += temp_score[index];
					temp_score = getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)	
					for index in xrange(int(numPoints)):		
						minus_scores[index] += temp_score[index];
	return (plus_scores, minus_scores);
	

def getTESProfile(coords, upstream_length, downstream_length, resolution, window_size, bed_vals):
	"""
	Draw profiles around a set of TESes 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary of TUCSC-type objects keyed by chrom
	bed-vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is pure count, not normalized at all
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		exit (1)
	numPoints = length/resolution;
	plus_scores = [];
	minus_scores = [];
    	for i in range(int(numPoints)):
        	plus_scores.append(0);
        	minus_scores.append(0);
	
	for chrom in coords.keys():
        	genes = coords[chrom];
        	if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = breakUpStrands(bed_vals[chrom]);
			
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			for g in genes:
				if plus.match(g.strand):
					temp_score = getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)
					for index in xrange(int(numPoints)):
						plus_scores[index] += temp_score[index];
					temp_score = getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)	
					for index in xrange(int(numPoints)):
						minus_scores[index] += temp_score[index];	
				elif minus.match(g.strand):
					## notice, swapped strands for tags because we're working on the crick strand now
					# ie, for genes on the minus strand, plus reads --> minus reads, minus-reads --> plus reads
					temp_score = getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)
					for index in xrange(int(numPoints)):
						plus_scores[index] += temp_score[index];
					temp_score = getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)	
					for index in xrange(int(numPoints)):		
						minus_scores[index] += temp_score[index];
	return (plus_scores, minus_scores);

def getTFBSProfile(coords, upstream_length, downstream_length, resolution, window_size, bed_vals):
	"""
	Draw profiles around a set of TFBSs. They do not have any orientations (I suppose). 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary keyed by chrom, with value being a list of locations.
	bed-vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is pure count, not normalized at all
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		sys.exit (1)
	numPoints = length/resolution;
	plus_scores = [];
	minus_scores = [];
    	for i in range(int(numPoints)):
        	plus_scores.append(0);
        	minus_scores.append(0);
	
	for chrom in coords.keys():
        	sites = coords[chrom];
        	if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = breakUpStrands(bed_vals[chrom]);
			
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			for site in sites:
				temp_score = getProfileNearPosition(site, '+', upstream_length, downstream_length, resolution, window_size, plus_starts)
				for index in xrange(int(numPoints)):
					plus_scores[index] += temp_score[index];
				temp_score = getProfileNearPosition(site, '+', upstream_length, downstream_length, resolution, window_size, minus_starts)	
				for index in xrange(int(numPoints)):
					minus_scores[index] += temp_score[index];
	return (plus_scores, minus_scores);


def output(upstream_length, resolution, plus_scores, minus_scores, normalization, outfilename):
	numPoints = len(plus_scores);
	assert (numPoints == len(minus_scores));
	outfile = open(outfilename, 'w');
	data_start = -1.0 * upstream_length
	xValues=[0.0]*numPoints;
	half_partition = int(resolution/2.0);
	for x in xrange(numPoints):
		xValues[x] = data_start + half_partition + x * resolution
		plus_norm_score = float(plus_scores[x]) / float(normalization);
		minus_norm_score = float(minus_scores[x]) / float(normalization);
		outline = str(xValues[x]) + "\t" + str(plus_norm_score) + "\t" + \
			str(minus_norm_score) + "\n";
		outfile.write(outline);
	outfile.close();
	return xValues;

def main(argv):
	parser = OptionParser()
	parser.add_option("-k", "--known_genes_file", action="store", type="string",
			dest="known_file", help="file with known genes", metavar="<file>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bedfile", help="file with tags in bed format", metavar="<file>")
	parser.add_option("-c", "--TypeOfSites", action="store", type="string",
			dest="type", help="TSS, TES, TFBS", metavar="<str>")
	parser.add_option("-o", "--outfile", action="store", type="string",
			dest="outfile", help="outfile name", metavar="<file>")
	parser.add_option("-n", "--normalization", action="store", type="float", dest="norm",
			help="additional normalization in addition to number of sites, number of reads per million and window_size per 1K")
	parser.add_option("-s", "--species", action="store", type="string",
			dest="species", help="species", metavar="<str>")
	parser.add_option("-u", "--UpstreamExtension", action="store", type="int",
			dest="upstreamExtension", help="UpstreamExtension", metavar="<int>")
	parser.add_option("-d", "--DownstreamExtension", action="store", type="int",
			dest="downstreamExtension", help="DownstreamExtension", metavar="<int>")
	parser.add_option("-r", "--resolution", action="store", type="int",
			dest="resolution", help="resolution of the profile, eg, 5", metavar="<int>")
	parser.add_option("-w", "--WindowSize", action="store", type="int",
			dest="window_size", help="window size for averaging. When window size > resolution, there is smoothing", metavar="<int>")			
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 20:
		parser.print_help()
		sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	#t0 = time.time()
	
	SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1')
	num_genes = 0
	num_tags = 0
	profiles = {}
	numPoints = float(opt.upstreamExtension + opt.downstreamExtension)/float(opt.resolution)
	print "Upstream extension: ", opt.upstreamExtension;
	print "Downstream extension: ", opt.downstreamExtension;
	print "Resolution:", opt.resolution;
	print "Scanning window size: ", opt.window_size;
	print "Number of Points", numPoints
	plus_score_profile = [0]*int(numPoints)
	minus_score_profile = [0]*int(numPoints)
	
	if (opt.type == "TSS"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			mycoords={}
			chrombed = chrom + '.bed1';
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()):
					if (len(coords[chrom]) > 0):
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = getTSSProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, bed_vals);		
	elif (opt.type == "TES"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			mycoords={}
			chrombed = chrom + '.bed1';
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()):
					if (len(coords[chrom]) > 0):
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = getTESProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, bed_vals);		
	elif (opt.type == "TFBS"):
		# Build coords
		# Here we are assuming that the file has the format chrom location + .....for each line
		# chrom is sline[0], location is sline[1]
		coords={};
		if(opt.known_file):	
			infile = open(opt.known_file, 'r');
			for line in infile:
				""" check to make sure not a header line """
				if not re.match("track", line):
					line = line.strip();
					sline = line.split();
					if sline[0] not in coords.keys():
						coords[sline[0]] = [];
					coords[sline[0]].append(atoi(sline[1]))
			infile.close();
		
		for chrom in chroms:
			mycoords={}
			chrombed = chrom + '.bed1';
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if chrom in coords.keys():
					if len(coords[chrom]) > 0:
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = getTFBSProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, bed_vals);		
	else:
		print "Only three types of locations are allowed: TSS, TES, TFBS"
		sys.exit(1); 
	
	for chrom in profiles.keys():
		(plus_scores, minus_scores) = profiles[chrom]
		assert ( int(numPoints) == len (plus_scores) )
		assert ( int(numPoints) == len (minus_scores) )
		for i in xrange(int(numPoints)):
			plus_score_profile[i] += plus_scores[i]
			minus_score_profile[i] += minus_scores[i]
	
	SeparateByChrom.cleanup(chroms, '.bed1')
        normalization = num_tags/1000000.0;
        normalization *= num_genes;
        normalization *= opt.window_size/1000.0;
        normalization *= opt.norm;
	print "Number of locations: ", num_genes;
	print "Number of reads: ", num_tags
	print "normalization = ", normalization
	
	xValues = output(opt.upstreamExtension, opt.resolution, plus_score_profile, minus_score_profile, normalization, opt.outfile);

	#print time.time() - t0, "seconds"
	
if __name__ == "__main__":
	main(sys.argv)


