#!/usr/bin/env python
"""
This module does a global rescaling of all the numerical columns, used for profiles, etc
"""

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
import Utility_extended

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="infile", metavar="<file>", help="name for input file, which is a multicolumn text file")
	parser.add_option("-c", "--column", action="store", type="int", dest="c", metavar="<int>", help="the index of the id column, 0-based")
	parser.add_option("-r", "--rescale_factor", action="store", type="float", dest="rescale_factor", metavar="<float>", help="the rescale factor that will be multiplied to the numbers in columns except c")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="name for output file")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	myrescale_factor = atof(opt.rescale_factor)
	
	infile = open(opt.infile,'r')
	outfile = open(opt.outfile, 'w')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			newline = []
			if (len(sline)>0):
				for index in xrange(len(sline)):
					if index != opt.c:
						new_value = (atof(sline[index]) * myrescale_factor);
						newline.append(str(new_value))
					else:
						newline.append(sline[index])
				outfile.write('\t'.join(newline)+'\n')
		else:
			outfile.write(line)
	infile.close()
	outfile.close()
	
if __name__ == "__main__":
	main(sys.argv)