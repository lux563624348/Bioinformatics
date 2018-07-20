#!/usr/bin/env python

"""
Assign islands to repetitive elements and Calculate the number of islands that are overlapping to REs
"""


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
from operator import itemgetter

try:
   import cPickle as pickle
except:
   import pickle

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')


sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements')

import UCSC_revised
import GenomeData
import SeparateByChrom
import Utility_extended
import RepElements
import get_read_count_on_REs

plus = re.compile("\+");
minus = re.compile("\-");


def assign_islands_to_REs(re_file_dir, re_file_name, chrom, chrom_length, island_list, upstream, downstream, min_re_length=0):
	"""
	islands are non-overlapping [start, end]
	returns {(start,end): [id]}
	
	"""
	currentdir = os.getcwd()
	os.chdir(re_file_dir)
	#{id:RepElement}
	#
	known_repelements = RepElements.KnownRepElements.initiate_from_file([chrom], re_file_name)
	#print len(known_repelements.rep_elements.keys())
	island_flags = [0 for island in island_list]
	for myid in known_repelements.rep_elements.keys():
		myelement = known_repelements.rep_elements[myid]
		#No matter whether positive or negative, genoStart < genoEnd
		if plus.match(myelement.strand):
			start = max(myelement.genoStart - upstream, 0)
			end = min(myelement.genoEnd + upstream, chrom_length)
		elif minus.match(myelement.strand):
			start = max(myelement.genoStart - downstream, 0)
			end = min(myelement.genoEnd + upstream, chrom_length)
		else:
			print myelement
			print "strand not recognized"
			exit(1)
		
		region = (start, end)
		(start_index, end_index) = Utility_extended.find_islands_overlapping_with_region(region, island_list) #returns [island]
		for index in range(start_index, end_index):
			island_flags [index] = 1 # These islands overlap with REs
		#print re_file_name, region, start_index, end_index
	os.chdir(currentdir)
	return island_flags
	
	
	
def main(argv):
	parser = OptionParser()
	
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island bed file")
	parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="file with RE tree in pickle format")
	parser.add_option("-l", "--RE_annotation_file_location", action="store", type="string", dest="RE_file_location", metavar="<file>", help="location of RE files named in repClass_repFamily_repName.txt")
	parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension", help="upstream extension from start", metavar="<int>")
	parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension", help="downstream extension from end",  metavar="<int>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	#Separate_by_chrom on bedfile
	lib_name = (opt.bedfile).split('/')[-1] # remove directory
	suffix = lib_name.split('.')[-1] # txt
	lib_name = lib_name.split('.')[0] 
	extension = "-" + lib_name +'.' + suffix +"1"
	if Utility_extended.fileExists(opt.bedfile):
		if Utility_extended.chrom_files_exist(chroms, extension) != 1:
			SeparateByChrom.separateByChrom(chroms, opt.bedfile,  extension)
	else:
		print bedfile, " is not found";
		sys.exit(1)
	
	
	
	print "\nLoad the RE tree to get the RE file names"
	re_tree = pickle.load(open(opt.RE_Tree, 'rb'))
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	
	total_num_islands = 0
	total_num_RE_islands = 0
	
	#cycle through chrom
	for chrom in chroms:
		# Get the islands
		island_list = []
		print chrom
		chrom_length = chrom_lengths[chrom]
		chrombed = chrom + extension
		if Utility_extended.fileExists(chrombed):		
			# load in each island
			inf = open(chrombed,'r')
			for line in inf:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					start = int(sline[1])
					end = int(sline[2])
					island_list.append( (start, end) )
			inf.close()	
			if Utility_extended.is_tuplelist_sorted(island_list, 0) != 1:
				island_list.sort(key = itemgetter[0]) # sort by start, assume non-overlapping
		else:
			print "%s can not be found" %chrombed
			
		island_flags = [0 for island in island_list]
		
		min_re_length = 10
		for reClass in re_tree.keys():
			for reFamily in re_tree[reClass].keys():
				for reName in re_tree[reClass][reFamily]:
					re_file_name = "_".join([reClass, reFamily, reName]) + ".txt"
					#print re_file_name
					this_island_flags = assign_islands_to_REs(opt.RE_file_location, re_file_name, chrom, chrom_length, island_list,  opt.upstream_extension, opt.downstream_extension, min_re_length) 
					#Collect the results into island_flags
					for i in xrange(len(this_island_flags)):
						if this_island_flags[i] == 1:
							island_flags[i] = 1
	
		print "There are %d island on %s" %(len(island_list),chrom)
		print "There are %d RE islands" %(sum(island_flags))
		total_num_islands += len(island_list)
		total_num_RE_islands += sum(island_flags)
	
	SeparateByChrom.cleanup(chroms, extension)
	print "There are %d islands" %(total_num_islands)
	print "There are %d RE islands" %(total_num_RE_islands)


	
if __name__ == "__main__":
	main(sys.argv)