#!/usr/bin/env python
"""
The islands can be overlapping
#
# 
The goal is to output 1) the read count on each island and 2) the island filtered reads

"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
from operator import itemgetter
import operator
import bisect

## get BED module
import BED_revised
from GenomeData import *
import SeparateByChrom
import Utility_extended

plus = re.compile("\+");
minus = re.compile("\-");

def tag_position(sline, shift):
	if plus.match(sline[5]):
		return atoi(sline[1]) + shift
	elif minus.match(sline[5]):
		return atoi(sline[2]) - 1 - shift


	
def find_reads_on_regions(read_file, regions, shift, outfile_name, boundary_extension = 0):
	"""
	regions: [BED3]. The regions can overlap
	read_file contains reads from only one chrom
	Return:
	read count on each island: {(start, end): rc}
	"""
	regions_rc = [] #[(region, rc)]
	tag_list = [] #[(position, line)]
	tag_flag_list = [] # Flags whether a read pass filtering.
	
	if Utility.fileExists(read_file):
		f = open(read_file,'r')
		for line in f:
			if not re.match("#", line):
				line = line.strip()
				sline = line.split()
				position = tag_position(sline, shift)
				tag_list.append((position, line))
		f.close()
	# sort according to position
	if (Utility_extended.is_listT_sorted(tag_list)==0):
		tag_list.sort(key=itemgetter(0))
	positions = [tag[0] for tag in tag_list]
	tag_flag_list = [0 for tag in tag_list]
	
	for region in regions:
		start = max(region.start - boundary_extension, 0)
		end = region.end + boundary_extension 
		assert (start<=end)
		start_ind = bisect.bisect_left(positions, start)
		end_ind = bisect.bisect_right(positions, end)
		for index in range(start_ind, end_ind): # These reads are on regions
			tag_flag_list[index] = 1
		rc = end_ind - start_ind
		regions_rc.append((region, rc))

	o = open(outfile_name, 'w')
	for i in xrange(len(tag_list)):
		if tag_flag_list[i] == 1:
			o.write( tag_list[i][1] + '\n') # tag_list[i][1] = line
	o.close();

	return regions_rc


def output_island_with_rc(island_rc, outfile_name):
	"""
	island_rc: [(island, rc)], island is a BED_annotated 3 object
	"""
	outf = open(outfile_name, "w")
	for item in island_rc:
		island = item[0]
		rc = item[1]
		outline = island.chrom + "\t" + island.start  + "\t" +  island.end + "\t" + str(rc) + "\t" + island.annotation + "\t"
		outf.write(outline)
	outf.close()
		
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store",
			  type="string", dest="species",
			  help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-a", "--rawbedfile", action="store",
			  type="string", dest="bedfile",
			  metavar="<file>",
			  help="raw data file in bed format")
	parser.add_option("-i", "--shift", action="store",
			  type="int", dest="shift",
			  metavar="<int>",
			  help="shift for finding the center of DNA fragment represented by the read")
	parser.add_option("-b", "--islandfile", action="store", type="string",
			  dest="islandbedfile", metavar="<file>",
			  help="island file")
	parser.add_option("-o", "--outfile", action="store", type="string",
			  dest="out_file", metavar="<file>",
			  help="filtered raw bed file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	islands = BED_revised.BED_annotated(opt.species, opt.islandbedfile, "BED3", 0);
	
	libName = (opt.bedfile).split('/')[-1] #remove directories
	libName = libName.split('.')[0] #remove .bed
	extension = "-" + libName +'.bed1'
	SeparateByChrom.separateByChrom(chroms, opt.bedfile, extension)
	for chrom in chroms:
		if chrom in islands.keys():
			outfile = chrom + "-filtered" + extension 
			# [(island, rc)]
			island_rc = filter_tags_by_islands(chrom + extension, islands[chrom], opt.shift, outfile, boundary_extension=0)
			output_island_with_rc(island_rc, chrom + "-islands" + extension)
			
	SeparateByChrom.combineAllGraphFiles(chroms, "-islands" + extension, "islands" + extension);
	SeparateByChrom.combineAllGraphFiles(chroms, "-filtered" + extension, opt.out_file);
	
	SeparateByChrom.cleanup(chroms, extension);
	SeparateByChrom.cleanup(chroms, "-filtered" + extension);
	SeparateByChrom.cleanup(chroms, "-islands" + extension);


if __name__ == "__main__":
	main(sys.argv)
