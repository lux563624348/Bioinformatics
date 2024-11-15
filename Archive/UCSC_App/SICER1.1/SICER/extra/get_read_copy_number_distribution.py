#!/usr/bin/env python
"""
Calculate the normalized read redundancy distribution given a bed file
"""
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import copy

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import GenomeData
import Utility_extended
import SeparateByChrom

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";
plus = re.compile("\+");
minus = re.compile("\-");



def count_redundancy_1chrom_single_strand_sorted(infile):
	'''
	infile must have bed format, can only contain reads from one chromosome and only one kind of strands (+/-). file must be pre-sorted by column2, then by column3.
	'''
	copy_number_histogram = {}
	
	current_start = 0
	current_end = 0
	current_count = 1 
	total = 0
	f = open(infile,'r')
	for line in f:
		if not re.match("#", line):
			total += 1
			line = line.strip()
			sline = line.split()
			start = atoi(sline[1])
			end = atoi(sline[2])
			if total == 1: #First read 
				current_start = start
				current_end = end
				current_count = 1
			if total >= 2:
				if start != current_start or end != current_end:
					# read copies exhaused, register the copy number counting
					if current_count in copy_number_histogram.keys():
						copy_number_histogram[current_count] += 1
					else:
						copy_number_histogram[current_count] = 1
					current_start = start
					current_end = end
					current_count = 1
				else:
					current_count += 1
	f.close()
	# The final entry
	if total > 0:
		if current_count in copy_number_histogram.keys():
			copy_number_histogram[current_count] += 1
		else:
			copy_number_histogram[current_count] = 1
	#print total
	#print copy_number_histogram
	return copy_number_histogram, total

def break_strand_and_tally(chrom):
	'''
	infile can only contain reads from one chromosome
	'''
	infile = chrom + ".bed1"
	outfile = chrom + ".bed2"
	
	# get plus reads and sort according to start and then end
	try:
		if os.system('grep [[:space:]]+ %s | sort -g -k 2,3 > plus.bed1' % (infile)):
			raise
	except: sys.stderr.write("+ reads do not exist in " + str(infile) + "\n");
	plus_read_copy_number_histogram, plus_total = count_redundancy_1chrom_single_strand_sorted('plus.bed1')
	
	
	# get minus reads and sort according to start and then end
	try:
		if os.system('grep [[:space:]]- %s | sort -g -k 2,3 > minus.bed1' % (infile)):
			raise
	except: sys.stderr.write("- reads do not exist in " + str(infile) + "\n");
	minus_read_copy_number_histogram, minus_total = count_redundancy_1chrom_single_strand_sorted('minus.bed1')
	
	read_copy_number_histogram = combine(plus_read_copy_number_histogram, minus_read_copy_number_histogram)
	
	total = plus_total + minus_total
	os.system('rm plus.bed1')
	os.system('rm minus.bed1')
	
	return read_copy_number_histogram, total
	
def combine(dic1, dic2):
	combined = {}
	mykeys = list( set(dic1.keys()).union( set(dic2.keys())) )
	for mykey in mykeys:
		if mykey in dic1.keys():
			combined[mykey] = dic1[mykey]
			if mykey in dic2.keys():
				combined[mykey] += dic2[mykey]
		else:
			combined[mykey] = dic2[mykey]
	return combined
	
	
	
def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--raw_bed_file", action="store", type="string",
                      dest="bed_file", help="raw bed file", metavar="<file>")    
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
		
	SeparateByChrom.separateByChrom(chroms, opt.bed_file, '.bed1')
	
	total_copy_number_distribution = {}
	total_total = 0
	for chrom in chroms:
		if (Utility_extended.fileExists(chrom + ".bed1")):
			copy_number_distribution, total = break_strand_and_tally(chrom)
			total_copy_number_distribution = copy.copy(combine(copy_number_distribution, total_copy_number_distribution))
			total_total += total
	
	#get total read count 
	total= 0
	for mykey in total_copy_number_distribution.keys():
		total += mykey * total_copy_number_distribution[mykey]
	print "There are %f reads in %s" %(total, opt.bed_file)
	#print total_total
	
	#Normalize
	normalization = sum(total_copy_number_distribution.values())
	for mykey in total_copy_number_distribution.keys():
		total_copy_number_distribution[mykey] /=normalization *1.0
	
	#Output
	Utility_extended.output_dic(total_copy_number_distribution, opt.out_file, "key", False)
	SeparateByChrom.cleanup(chroms, '.bed1')
	
 
if __name__ == "__main__":
	main(sys.argv)
