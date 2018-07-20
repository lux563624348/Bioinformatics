#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


modifications_38 = ['H2AK5ac','H2AK9ac','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K91ac']

modifications = ['H2AK5ac','H2AK9ac','H2AK127me1','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K23me2','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K79me2','H4K91ac','H4R3me2']

acetylations = ['H2AK5ac','H2AK9ac','H2BK5ac','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K9ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K36ac','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K91ac']

methylations = ['H2AK127me1','H2BK5me1','H3K4me1','H3K4me2','H3K4me3','H3K9me1','H3K9me2','H3K9me3','H3K23me2','H3K27me1','H3K27me2','H3K27me3','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K20me1','H4K20me3','H4K79me2','H4R3me2']


def gene_modification_denary_set(datafile):
	'''returns a dictionary, whose keys are geneIDs, and values are the denary integers converted from the binary strings of the gene modification patterns'''
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene = sline[0]
			del sline[0]
			result[gene] = atoi(''.join(sline),2)
	file.close()
	return result


def average(List):
	total = 0.0
	for item in List:
		total += item
	return total/len(List)


def middle(List):
	if len(List)!=0:
		if len(List)%2 == 1:
			return List[len(List)/2]
		else:
			return (List[len(List)/2-1]+List[len(List)/2])/2.0
	else: return 0.0


def standard_deviation(List):
	total = 0.0
	squared_total = 0.0
	for item in List:
		total += item
		squared_total += pow(item,2)
	return sqrt(squared_total/len(List) - pow(total/len(List),2))


def merge_lists(List1, List2):
	List = List1
	for item in List2:
		if not item in List1:
			List.append(item)
	return List


def read_complete_binary_file(datafile):
	'''returns a dictionary, whose keys are geneIDs. the value is a list storing the sequence of binary numbers 0 or 1 of the absence/occurance of the modification, ordering corresponding to the modifications list'''
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


def get_modification_number_for_genes(datadic, genelist):
	'''return a dictionary, geneID as keys, value tells how many modifications this gene has among all 42 markers'''
	result = {}
	for gene in genelist:
		if datadic.has_key(gene):
			n = 0
			for digit in datadic[gene]:
				n += digit
			result[gene] = n
	return result


def get_complete_modification_level_list(filetail, path):
	'''There should be a set of files whose names start with modification and have a common tail. returns a list, elements in the order of modifications list, of a dictionary of gene modification levels'''
	List = []
	for mod in modifications_38:
		List.append(get_expression_data_dic(path+'/'+mod+filetail, 1))
	return List


def get_modification_gene_fraction_list(datadic):
	'''from patterns, find the fraction of genes that have each modifications'''
	List = [0.0]*len(datadic[datadic.keys()[0]])
	for gene in datadic.keys():
		for i in range(0, len(List)):
			List[i] += datadic[gene][i]
	for i in range(0, len(List)):
		List[i] = List[i] / float(len(datadic.keys()))
	return List


def get_modification_gene_number_list(datadic):
	'''from patterns, find the fraction of genes that have each modifications'''
	List = [0.0]*len(datadic[datadic.keys()[0]])
	for gene in datadic.keys():
		for i in range(0, len(List)):
			List[i] += datadic[gene][i]
	return List


def get_gene_list_with_modification(datadic, genelist, number, value):
	'''from a list of genes, with their modification pattern data dictionary, returns a list of geneIDs, whose modification of number is value(0 or 1)'''
	List = []
	for gene in genelist:
		if datadic.has_key(gene):
			if datadic[gene][number] == value:
				List.append(gene)
	return List


def get_gene_list_with_modification_from_all(datadic, number, value):
	'''from modification pattern data dictionary, returns a list of geneIDs, whose modification of number is value(0 or 1)'''
	return get_gene_list_with_modification(datadic, datadic.keys(), number, value)


def get_gene_total_modification_level_sum_Log(modificationlevellist, genelist):
	'''returns a dictionary, whose keys are geneIDs from genelist, the value is the sum of Log10 value of all non-zero modification levels of the gene'''
	dic = {}
	for gene in genelist:
		total = 0.0
		for modification in modificationlevellist:
			assert modification.has_key(gene)
			if modification[gene] > 0.1:
				total += log(modification[gene],10)
		dic[gene] = total
	return dic


def get_data_binary_dic_subset(datadic, List):
	'''returns the subdictionary of datadic whose keys are in List'''
	dic = {}
	for item in List:
		if datadic.has_key(item):
			dic[item] = datadic[item]
	return dic


def get_expression_data_dic(expressionfile, colum):
	'''returns a dictionary with geneIDs as keys, expression value as values. colum is the number of colum -1 where the expression data are in the file'''
	file = open(expressionfile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			dic[sline[0]] = atof(sline[colum])
	file.close()
	return dic


def get_expression_log_data_dic(expressionfile, colum):
	'''similar to get_expression_data_dic, but returns the log10 of expression value'''
	file = open(expressionfile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			dic[sline[0]] = log(atof(sline[colum]),10)
	file.close()
	return dic

def get_positive_log_data_dic(expressionfile, colum):
	'''similar to get_expression_data_dic, but returns the log10 of expression value, 0 if number < 1'''
	file = open(expressionfile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			n = atof(sline[colum])
			if n < 1.0:
				dic[sline[0]] = 0.0
			else:
				dic[sline[0]] = log(n,10)
	file.close()
	return dic


def get_expression_list_from_gene_list(dic, List):
	'''returns the list of expression values of the genes in the input list'''
	datalist = []
	for gene in List:
		if dic.has_key(gene):
			datalist.append(dic[gene])
	return datalist


def get_pattern_gene_set(datadic, patternlist):
	'''returns a list of the genes who have the modification pattern as described in patternlist'''
	genelist = datadic.keys()
	for i in range(0, len(patternlist)):
		if patternlist[i] == 0 or patternlist[i] == 1:
			resultlist = get_gene_list_with_modification(datadic, genelist, i, patternlist[i])
			genelist = resultlist
	return genelist


def get_pattern_list_from_file(patternfile, length):
	'''input file should be a one line string of numbers seperated by space or tab, with '1' (occurance), '0' (absence), or '?' (dont care whether this modification exists or not). Returns a 42 long list with 1(occurance), 0(absence), or -1(dont care)'''
	pattern_list = [-1]*length
	file = open(patternfile, 'r')
	for line in file:
		line = line.strip()
		sline = line.split()
		assert len(sline) == length
		for i in range(0, len(sline)):
			if not re.match("?", sline[i]):
				pattern_list[i] = atoi(sline[i])
	file.close()
	return pattern_list


def get_pattern_list_from_string(string, length):
	'''to convert a binary tring to a list of binary number of 0/1'''
	pattern_list = [-1]*length
	assert len(string) == length
	for i in range(0, len(string)):
		if string[i] == '1' or string[i] == '0':
			pattern_list[i] = atoi(string[i])
	return pattern_list


def average_expression_of_genes(datadic, genelist):
	'''returns the mean value of the expression levels of the genes in genelist'''
	return average(get_expression_list_from_gene_list(datadic, genelist))


def expression_deviation_of_genes(datadic, genelist):
	return standard_deviation(get_expression_list_from_gene_list(datadic, genelist))


def get_mean_expression_change_with_modification(datadic, expressiondic, modification):
	'''returns the difference of the mean expression level of the genes with modification and that of the genes without modification'''
	if modification in modifications:
		index = modifications.index(modification)
		genelist1 = get_gene_list_with_modification_from_all(datadic, index, 1)
		genelist0 = get_gene_list_with_modification_from_all(datadic, index, 0)
		average1 = average_expression_of_genes(expressiondic, genelist1)
		average0 = average_expression_of_genes(expressiondic, genelist0)
		return average1 - average0
	else:
		print "Modification input Error"


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--completedata", action="store", type="string", dest="infile", metavar="<file>", help="input whole data set file name")
	parser.add_option("-k", "--expressionfile", action="store", type="string", dest="expressionfile", metavar="<file>", help="input gene expression file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	datadic = read_complete_binary_file(opt.infile)
	expressiondic = get_expression_log_data_dic(opt.expressionfile, 12)
	f = open(opt.output,'w')
	for mod in modifications:
		f.write(mod + '\t' + str(get_mean_expression_change_with_modification(datadic, expressiondic, mod)) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)