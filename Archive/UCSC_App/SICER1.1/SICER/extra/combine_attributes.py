#!/usr/bin/env python

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import Utility

def main(argv):
	# Output
	# Gene Name    read count    read count normalized by total per million
	parser = OptionParser()
	parser.add_option("-a", "--fA", action="store", type="string", dest="fA", metavar="<file>", help="file with attributes set A")
	parser.add_option("-b", "--fB", action="store", type="string", dest="fB", metavar="<file>", help="file with attributes set B")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for ids with combined attributes")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
	if not Utility.fileExists(opt.fA):
		print opt.fA, " not found";
		sys.exit(1)
	if not Utility.fileExists(opt.fB):
		print opt.fB, " not found";
		sys.exit(1)
		
	
	genes = {}; # a BED object of BED_GRAPH elements
	

	f = open(opt.fA,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > 1: 
				genes[sline[0]] = sline[1:];
	f.close();
	
	f = open(opt.fB,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (sline[0] in genes.keys()) and len(sline) > 1:
				genes[sline[0]].extend(sline[1:]);
	f.close();			
				
	f= open(opt.out_file, 'w')
	for item in genes.keys():
		sline = "\t".join(genes[item]);
		line = item + "\t" + sline + "\n";
		f.write(line)
	f.close()
	
	print "Total number of ",  opt.region_type , ": ", len(genes.keys());


if __name__ == "__main__":
	main(sys.argv)