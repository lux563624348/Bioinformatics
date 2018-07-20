#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import copy

import gene_set_manipulation

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();



def read_complete_binary_file(datafile):
	'''
	returns a dictionary, whose keys are geneIDs. the value is a list storing the sequence of binary numbers 0 or 1 of the absence/occurance of the feature, ordering corresponding to the feature list
	'''
	file = open(datafile,'r')
	result = {}
	for line in file:
		List = []
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene = sline[0]
			for i in range(1, len(sline)):
				List.append(atoi(sline[i]))
			result[gene] = List
	file.close()
	return result

def get_pattern_gene_set(data_dic, patternlist):
	'''
	patternlist = [0 or 1 or N]
	returns a list of the genes who have the feature pattern as described in patternlist
	'''
	genelist = data_dic.keys()
	for i in range(0, len(patternlist)):
		if patternlist[i] == 0 or patternlist[i] == 1:
			resultlist = get_gene_list_with_feature(data_dic, genelist, i, patternlist[i])
			genelist = copy.copy(resultlist)
	return genelist

def get_gene_list_with_feature(data_dic, genelist, index, value):
	'''
	data_dic = {name:[feature values]}
	value:
		0: absent
		1: present
	from a list of genes, with their feature pattern data dictionary, returns a list of geneIDs, whose feature of index is value(0 or 1)
	'''
	List = []
	for gene in genelist:
		if data_dic.has_key(gene):
			if data_dic[gene][index] == value:
				List.append(gene)
	return List

def get_pattern_list_from_string(string):
	'''
	to convert a string of 0/1/N to a list of binary number of 0/1
	pattern_list = [0 or 1 or -1], for those places with N it is a -1 
	'''
	length = len(string)
	pattern_list = [-1]*length
	assert len(string) == length
	for i in range(0, len(string)):
		if string[i] == '1' or string[i] == '0':
			pattern_list[i] = atoi(string[i])
	return pattern_list
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--completedata", action="store", type="string", dest="infile", metavar="<file>", help="input digitized data file name")
	parser.add_option("-k", "--annotation_file", action="store", type="string", dest="annotation_file", metavar="<file>", help="input gene expression file name")
	parser.add_option("-p", "--pattern", action="store", type="string", dest="pattern", metavar="<str>", help="combinatorial patterning string, eg, 0011NN")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	data_dic = read_complete_binary_file(opt.infile) #{name:[feature values]}
	pattern_list = get_pattern_list_from_string(opt.pattern)
	genelist = get_pattern_gene_set(data_dic, pattern_list)
	print "There are %d genes with the required pattern" %len(genelist)
	gene_set_manipulation.output_UCSCsubset_in_file (opt.annotation_file, genelist, opt.output)


if __name__ == "__main__":
	main(sys.argv)