#!/usr/bin/env python
"""
Digitize read count value, used in conjuction with ... for pattern analysis

"""
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import gene_set_manipulation

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def get_values(datafile):
	"""
	result:{name:rc}
	"""
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result[sline[0]] = atof(sline[1])
	file.close()
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputFileNameSuffix", action="store", type="string", dest="infile_suffix", metavar="<STRING>", help="input file name suffix")
	parser.add_option("-l", "--featurelist", action="store", type="string", dest="feature_list_file", metavar="<file>", help="list of all feature in a file")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	feature_list = gene_set_manipulation.get_gene_list(opt.feature_list_file, 0)
	
	complete_data = {}
	for feature in feature_list:
		complete_data[feature] = get_values(feature + opt.infile_suffix)
	
	f = open(opt.output,'w')
	f.write('#Gene Name')
	for feature in feature_list:
		f.write('\t' + feature)
	for gene in complete_data[feature_list[0]].keys():
		f.write('\n'+gene)
		for feature in feature_list:
			if complete_data[feature][gene] > 0.000001:
				f.write('\t1')
			else:
				f.write('\t0')
	f.close()


if __name__ == "__main__":
	main(sys.argv)