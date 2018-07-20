#!/usr/bin/env python

# Calculate read counts on chroms and normalize the rc


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/JunZhu/modules')

import get_total_tag_counts
import get_chrom_length	

def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--fastafile", action="store", type="string", dest="fasta_file", metavar="<file>", help="fasta file for the sequences")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedFile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	
	chrom_lengths = get_chrom_length.get_chrom_lengths(opt.fasta_file)
	print "There are %d RE species" %len(chrom_lengths)
	
	chrom_rc = {}
	
	f = open(opt.out_file, 'w')
	outline = "# Name"  + "\t" + "RC" + "\t" + "RPM" + "\t"  + "RPKM" + '\n'
	f.write(outline)
	
	inf = open(opt.bedFile,"r")
	for line in inf:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			chrom = sline[0]
			if chrom in chrom_rc.keys():
				chrom_rc[chrom] += 1.0
			else:
				chrom_rc[chrom] = 1.0
	inf.close()
	
	#chroms with reads
	chroms = list(set(chrom_rc.keys()) & set(chrom_lengths.keys()))
	
	total_length = sum(chrom_lengths.values())
	totalcount = get_total_tag_counts.get_total_tag_counts(opt.bedFile) * 1.0
	#print totalcount
	
	basel_rpkm = (1000000.0)/(total_length/1000.0)
	
	print "Basel RPKM for chroms with reads is %f" %basel_rpkm 
	
	for chrom in chroms:
		rc = chrom_rc[chrom]
		rpm = (rc/totalcount)*1000000
		rpkm = rpm * 1000.0 / chrom_lengths[chrom]
		outline = chrom + "\t" + str(rc) + "\t" + str(rpm) + "\t" + str(rpkm) +"\n"
		f.write(outline)
	f.close()

if __name__ == "__main__":
	main(sys.argv)