#!/usr/bin/env python

"""
03/2011

This is a template for the analysis of tag distribution with respect
to a set of points, such as the TSSs of known genes.

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
import GenerateProfileAroundLocations


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	desc="""This is a template for the analysis of aggretated tag distribution with respect to a set of points, such as the TSSs of known genes, with the reads from the two strands combined."""
	parser = OptionParser(description=desc)
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
	parser.add_option("-i", "--shift", action="store", type="int", dest="shift", help="amount of shift for each read to get the center of the fragment ", metavar="<int>")		
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 22:
		parser.print_help()
		sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species];
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
	profiles = {}
	numPoints = float(opt.upstreamExtension + opt.downstreamExtension)/float(opt.resolution)
	print "Upstream extension: ", opt.upstreamExtension;
	print "Downstream extension: ", opt.downstreamExtension;
	print "Resolution:", opt.resolution;
	print "Scanning window size: ", opt.window_size;
	print "Number of Points", numPoints
	score_profile = [0]*int(numPoints)
	
	if (opt.type == "TSS"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			mycoords={}
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()):
					if (len(coords[chrom]) > 0):
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = GenerateProfileAroundLocations.getTSSPMProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.shift, opt.shift, bed_vals)		
	elif (opt.type == "TES"):
		coords = UCSC.KnownGenes(opt.known_file);
		for chrom in chroms:
			mycoords={}
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if (chrom in coords.keys()):
					if (len(coords[chrom]) > 0):
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = GenerateProfileAroundLocations.getTESPMProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.shift, opt.shift, bed_vals)	
							
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
			chrombed = chrom + extension;
			if Utility.fileExists(chrombed):
				bed_vals = {};
				bed_vals = BED.BED(opt.species, chrombed, "BED2");
				num_tags += bed_vals.getNumVals()
				if chrom in coords.keys():
					if len(coords[chrom]) > 0:
						num_genes += len(coords[chrom])
						mycoords[chrom] = coords[chrom];
						profiles[chrom] = GenerateProfileAroundLocations.getTFBSPMProfile(mycoords, opt.upstreamExtension, opt.downstreamExtension, opt.resolution, opt.window_size, opt.shift, opt.shift, bed_vals)		
	else:
		print "Only three types of locations are allowed: TSS, TES, TFBS"
		sys.exit(1); 
	
	for chrom in profiles.keys():
		assert ( int(numPoints) == len (profiles[chrom]) )
		for i in xrange(int(numPoints)):
			score_profile[i] += profiles[chrom][i]
	
	SeparateByChrom.cleanup(chroms, extension)
	normalization = num_tags/1000000.0;
	normalization *= num_genes;
	normalization *= opt.window_size/1000.0;
	normalization *= opt.norm;
	print "Number of locations: ", num_genes;
	print "Number of reads: ", num_tags
	print "normalization = ", normalization
	
	GenerateProfileAroundLocations.outputPMP(opt.upstreamExtension, opt.resolution, score_profile, normalization, opt.outfile)


	#print time.time() - t0, "seconds"
	
if __name__ == "__main__":
	main(sys.argv)


