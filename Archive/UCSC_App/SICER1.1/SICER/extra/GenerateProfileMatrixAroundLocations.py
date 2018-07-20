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
import pylab

import Utility
import BED
import UCSC
import GenomeData 
import SeparateByChrom
import GenerateProfileAroundLocations


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def getTSSPMProfileMatrix(coords, upstream_length, downstream_length, resolution, window_size, pshift, mshift, bed_vals):
	"""
	Draw profiles around a set of TSSes 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary of UCSC-type objects keyed by chrom
	bed_vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is a list of profiles, each item in the list represents a profile around one TSS.  Pure count, not normalized at all
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		exit (1)
	numPoints = length/resolution;

	scores = [0] * int(numPoints);
	scoredic = {}
	
	for chrom in coords.keys():
		genes = coords[chrom];
		if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = GenerateProfileAroundLocations.breakUpStrands(bed_vals[chrom]);
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			plus_starts = [ item + pshift for item in plus_starts ]
			minus_starts = [ item - mshift for item in minus_starts ]
			
			for g in genes:
				if plus.match(g.strand):
					plus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)
					minus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)	
				elif minus.match(g.strand):
					## notice, swapped strands for tags because we're working on the crick strand now
					# ie, for genes on the minus strand, plus reads --> minus reads, minus-reads --> plus reads
					plus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)
					minus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)
					#print "\t".join(annotation)
					#print plus_scores
					#print minus_scores
					#print scoreMatrix
				scores = [sum(pair) for pair in zip(plus_scores, minus_scores)]
				# This causes error as the scores are kept as reference in scoreMatrix, mutation in scores cause changes to previous entries. Ways to resolve it include the above, or set scores =[] within the loop. When [] is encountered, a rebinding of scores with a new copy happens. 
				#for index in xrange(int(numPoints)):
				#	scores[index] = plus_scores[index] + minus_scores[index]
				#print scores
				#print scoreMatrix
				scoredic[g.name] = scores
				
	return scoredic;

def getTESPMProfileMatrix(coords, upstream_length, downstream_length, resolution, window_size, pshift, mshift, bed_vals):
	"""
	Draw profiles around a set of TESes 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary of TUCSC-type objects keyed by chrom
	bed_vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is a list of profiles, each item in the list represent a profile around one TES.  Pure count, not normalized at all
		  a list of corresponding annotations.
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		exit (1)
	numPoints = length/resolution;
	
	scores = [0] * int(numPoints); 
	scoredic = {}
	
	for chrom in coords.keys():
		genes = coords[chrom];
		if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = GenerateProfileAroundLocations.breakUpStrands(bed_vals[chrom]);
			
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			plus_starts = [ item + pshift for item in plus_starts ]
			minus_starts = [ item - mshift for item in minus_starts ]
			
			for g in genes:
				if plus.match(g.strand):
					annotation = [g.name, g.chrom, g.strand, str(g.txEnd)]
					plus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)
					minus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txEnd, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)	
				elif minus.match(g.strand):
					annotation = [g.name, g.chrom, g.strand, str(g.txStart)]
					## notice, swapped strands for tags because we're working on the crick strand now
					# ie, for genes on the minus strand, plus reads --> minus reads, minus-reads --> plus reads
					plus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, minus_starts)
					minus_scores = GenerateProfileAroundLocations.getProfileNearPosition(g.txStart, g.strand, upstream_length, downstream_length, resolution, window_size, plus_starts)	
				scores = [sum(pair) for pair in zip(plus_scores, minus_scores)]
				scoredic[g.name] = scores
	return scoredic;	

def getTFBSPMProfileMatrix(coords, upstream_length, downstream_length, resolution, window_size, pshift, mshift, bed_vals):
	"""
	Draw profiles around a set of TFBSs. They do not have any orientations (I suppose). 
	resolution is the partition size. 
		When window size is bigger than the resolution, there is a sliding window averaging.
		When window size = resolution, simple partitioning
	coord is a dictionary keyed by chrom, with value being a list of locations.
	bed_vals is a dictionary of reads keyed by chrom
		Tags also have orientation, they are taken into account elsewhere.
	output is a list of profiles, each item in the list represent a profile around one TES.  Pure count, not normalized at all
		  a list of corresponding annotations.
	
	"""
	length = upstream_length + downstream_length;
	if (length%resolution != 0):
		print "Please choose the resolution commensurate with the length of the regions"
		sys.exit (1)
	numPoints = length/resolution;
	scores = [0] * int(numPoints); 
	scoredic = {}	
	
	for chrom in coords.keys():
		sites = coords[chrom];
		if chrom in bed_vals.keys():
			(plus_starts, minus_starts) = GenerateProfileAroundLocations.breakUpStrands(bed_vals[chrom]);
			if not Utility.is_list_sorted(plus_starts):
				plus_starts.sort();
			if not Utility.is_list_sorted(minus_starts):
				minus_starts.sort();
			plus_starts = [ item + pshift for item in plus_starts ]
			minus_starts = [ item - mshift for item in minus_starts ]
			for site in sites:
				name = chrom + str(site)
				annotation = [name, chrom, '+', str(site)]
				plus_scores= GenerateProfileAroundLocations.getProfileNearPosition(site, '+', upstream_length, downstream_length, resolution, window_size, plus_starts)
				minus_scores= GenerateProfileAroundLocations.getProfileNearPosition(site, '+', upstream_length, downstream_length, resolution, window_size, minus_starts)
				scores = [sum(pair) for pair in zip(plus_scores, minus_scores)]	
				scoredic[g.name] = scores
	return scoredic;	

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
			help="additional normalization in addition to number of reads per million and window_size per 1K")
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
	parser.add_option("-p", "--plusReadShift", action="store", type="int",
			dest="pshift", help="plusReadShift", metavar="<int>")
	parser.add_option("-m", "--minusReadShift", action="store", type="int", 
			dest="mshift", help="minusReadShift", metavar="<int>")		
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 24:
		parser.print_help()
		sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	#t0 = time.time()
	libName = (opt.bedfile).split('/')[-1]
	libName = libName.split('.')[0]
	extension = "-" + libName +'.bed1'
	SeparateByChrom.separateByChrom(chroms, opt.bedfile, extension)
	
	num_genes = 0
	num_tags = 0
	numPoints = float(opt.upstreamExtension + opt.downstreamExtension)/float(opt.resolution)
	print "Upstream extension: ", opt.upstreamExtension;
	print "Downstream extension: ", opt.downstreamExtension;
	print "Resolution:", opt.resolution;
	print "Scanning window size: ", opt.window_size;
	print "Number of Points", numPoints
	
	all_genes_scores = {} #{name:[]}
	
	if (opt.type == "TSS"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				scoredic={}
				mycoords={}
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()) and (len(coords[chrom]) > 0):
					num_genes += len(coords[chrom])
					mycoords[chrom] = coords[chrom];
					scoredic = getTSSPMProfileMatrix(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.pshift, opt.mshift, bed_vals)
					all_genes_scores.update(scoredic)
					#print annotations
					#print scoreMatrix
					#print score_profiles
	elif (opt.type == "TES"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				scoredic={}
				mycoords={}
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()) and (len(coords[chrom]) > 0):
					num_genes += len(coords[chrom])
					mycoords[chrom] = coords[chrom];
					scoredic = getTESPMProfileMatrix(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.pshift, opt.mshift, bed_vals)
					all_genes_scores.update(scoredic)
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
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				scoredic={}
				mycoords={}
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if chrom in coords.keys() and len(coords[chrom]) > 0:
					num_genes += len(coords[chrom])
					mycoords[chrom] = coords[chrom];
					scoredic = getTFBSPMProfileMatrix(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.pshift, opt.mshift, bed_vals)
					all_genes_scores.update(scoredic)
	else:
		print "Only three types of locations are allowed: TSS, TES, TFBS"
		sys.exit(1); 
	
	SeparateByChrom.cleanup(chroms, extension)
	
	normalization = num_tags/1000000.0;
	normalization *= opt.window_size/1000.0;
	normalization *= opt.norm
	
	outFile = open(opt.outfile, 'w')
	# export the normalized result to a file. 
	for mykey in all_genes_scores.keys():
		outline = str(mykey) + "\t" + "\t".join( [str(item/normalization) for item in all_genes_scores[mykey]] ) + '\n'
		outFile.write(outline)
	outFile.close()
	
	
	print "Number of locations: ", num_genes;
	print "Number of reads: ", num_tags
	print "normalization = ", normalization
	
	# Testing
	overall_profile = [0] * int(numPoints);
	for mykey in all_genes_scores.keys():
		assert (len(all_genes_scores[mykey]) == int(numPoints))
		for j in xrange(int(numPoints)):
			overall_profile[j] += (all_genes_scores[mykey])[j]/normalization
	overall_profile = [item/float(num_genes) for item in overall_profile]
	#for item in overall_profile:
	#	print item;
	pylab.clf()
	pylab.plot(overall_profile,"b")
	pylab.savefig("Overall_profile.png", format='png')
	#print time.time() - t0, "seconds"
	
if __name__ == "__main__":
	main(sys.argv)


