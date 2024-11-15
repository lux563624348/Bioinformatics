#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter

## get BED module
import BED;
import UCSC_revised;
import Utility_extended


def get_gene_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	file = open(gene_file,'r')
	gene_list = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline)>0):
				gene_list.append(sline[c])
	file.close()
	return gene_list


def get_float_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[c]))
	file.close()
	return List


def get_int_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(int(atof(sline[c])))
	file.close()
	return List

			
def get_conversion_table(gene_file, origin, end):
	"""
	origin and end are the 0-based column numbers 
	Return a dictionary
	
	"""
	conversion_table = {};
	file = open(gene_file,'r')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > max(origin, end):
				conversion_table[sline[origin]] = sline[end];
	file.close()
	return conversion_table;


def get_gene_float_dic(expressionfile, colum):
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
	
	
	
def find_redundant_genes(gene_list):
	"""
	Return the list of names
	"""
	redundance_list=[];
	gene_list.sort();
	if len(gene_list)>=2:
		if gene_list[0] == gene_list[1]:
			redundance_list.append(gene_list[1]);
		for i in range(2, len(gene_list)):
			if gene_list[i] == gene_list[i-1] and gene_list[i-1] != gene_list[i-2]:
				redundance_list.append(gene_list[i]);
	return redundance_list;


def find_redundancy(gene_list):
	"""
	Return the list of (names, redundancy)
	"""
	redundancy_list=[]
	if len(gene_list) > 1:
		gene_list.sort()
	if len(gene_list) == 1:
		return [(gene_list[0], 1)]
	if len(gene_list) >= 2:
		current_gene = gene_list[0]
		redundancy = 1
		for i in range(1, len(gene_list)):
			if gene_list[i] == gene_list[i-1]:
				redundancy += 1
			else:
				redundancy_list.append((current_gene, redundancy))
				current_gene = gene_list[i]
				redundancy = 1
	return redundancy_list


def find_unique_genes (gene_list):
	unique_list=[];
	if gene_list !=[]:
		gene_list.sort();
		unique_list.append(gene_list[0]);
		for i in range(1, len(gene_list)):
			if gene_list[i] != gene_list[i-1]:
				unique_list.append(gene_list[i]);
	return unique_list;



def find_unique_only_list (gene_list):
	unique_list=[];
	if gene_list !=[]:
		gene_list.sort();
		if gene_list[1] != gene_list[0]:
			unique_list.append(gene_list[0]);
		for i in range(1, len(gene_list)-1):
			if gene_list[i] != gene_list[i-1] and gene_list[i] != gene_list[i+1]:
				unique_list.append(gene_list[i]);
		if gene_list[len(gene_list)-1] != gene_list[len(gene_list)-2]:
			unique_list.append(gene_list[len(gene_list)-1]);
	return unique_list;


def gene_comparison(List1, List2):
	"""
	This module is used to compare two sets of gene names and find the same and different genes between them. The input are two files, with the RefSeq IDs of the genes in the first column of the first file and in the second column of the second file. The output are three files of gene RefSeq lists, which are the same genes and the different genes for each one.
	""" 
	same = 'shared';
	diff1 = 'only in 1';
	diff2 = 'only in 2';
	
	result_lists = {}
	result_lists[same] = list(set(List1) & set(List2))
	result_lists[diff1] = list(set(List1) - set(List2))
	result_lists[diff2] = list(set(List2) - set(List1))
	return result_lists



def get_gene_dictionary(IDfile, c, c1 = 0):
	'''
	'''
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > max(c,c1):
				refseq[sline[c1]] = sline[c];
	file.close()
	return refseq


def convertRefSeqID(genes, conversion_table):
	"""
	RefSeq is the conversion table.
	"""
	for gene in genes:
		if conversion_table.has_key(gene['name']):
			gene['name'] = conversion_table[gene['name']]
	return genes


def convertID(list, conversion_table):
	"""
	conversion table is a dic
	"""
	for gene in list:
		if conversion_table.has_key(gene):
			gene = conversion_table[gene]
	return list


def convert_geneID_in_file(gene_file, geneID_column, conversion_table,  outfile):
	"""
	"""
	file = open(gene_file,'r')
	ofile = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_name = sline[geneID_column];
			#print gene_name;
			if conversion_table.has_key(gene_name):
				sline[geneID_column] = conversion_table[gene_name];	
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: ofile.write(line);
	ofile.close();
	file.close();


def output_UCSCsubset_in_file (known_gene_file, specific_gene_list, outfile):
	"""
	gene_file: is a file of tabular format, must be a UCSC file with name in the 0th column
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be extracted from gene_file and written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	"""
	output_subset_in_file (known_gene_file, 0, specific_gene_list, outfile);
			
	
def output_subset_in_file (gene_file, gene_name_column, specific_gene_list, outfile):
	"""
	gene_file: is a file of tabular format, could be a UCSC file or an expression file
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be extracted from gene_file and written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	This will not take care of redundancy in gene_file
	"""
	file = open(gene_file,'r');
	ofile = open(outfile,'w');
	print "There are ", len(specific_gene_list), " in your list";
	unique_gene_list = find_unique_genes(specific_gene_list);
	print len(unique_gene_list)," out of ", len(specific_gene_list), " genes in your list are unique";
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			if gene_name_column < len(sline):
				for item in unique_gene_list:
					#if re.match(sline[gene_name_column], item):
					if sline[gene_name_column] == item:
						outline = '\t'.join(sline) +'\n';
						ofile.write(outline);
						break; #Once found a match,no need to look further
	file.close();
	ofile.close();
	
def join(IDs, c1, file1, c2, file2, outfile):
	# IDs are shared by the two files.
	mylist_1 = build_list(IDs, c1, file1)
	mylist_2 = build_list(IDs, c2, file2)
	outlist = join_lists(mylist_1, mylist_2)
	ofile = open(outfile,'w')
	for item in outlist:
		out = [item[0]] + item[1]
		outline = "\t".join([str(i) for i in out]) + '\n'
		ofile.write(outline)
	ofile.close()
	
def join_lists(mylist_1, mylist_2):
	"""
	input: lists of (key, annotation)
	output: list of (key, annotation)
	
	"""
	if  Utility_extended.is_listT_sorted(mylist_1) != 1:
		mylist_1 = sorted(mylist_1, key=itemgetter(0))
	mylist_1 = remove_redundancy(mylist_1)
	#print len(mylist_1)
	if Utility_extended.is_listT_sorted(mylist_2) != 1:
		mylist_2 = sorted(mylist_2, key=itemgetter(0))
	mylist_2 = remove_redundancy(mylist_2)
	#print len(mylist_2)
	mylist_1_IDs = [i[0] for i in mylist_1]
	mylist_2_IDs = [i[0] for i in mylist_2]
	
	outlist = []
	if len(mylist_1) <= len(mylist_2):
		for item in mylist_1:
			ID = item[0]
			if ID in mylist_2_IDs:
				index = mylist_2_IDs.index(ID)
				out = item[1] + mylist_2[index][1]
				outlist.append((ID, out))
	else:
		for item in mylist_2:
			ID = item[0]
			if ID in mylist_1_IDs:
				index = mylist_1_IDs.index(ID)
				out = mylist_1[index][1] + item[1]
				outlist.append((ID, out))
	return 	outlist		
	
def join_dics(mydic_1, mydic_2):
	outdic = {}
	intersection = Utility_extended.intersection(mydic_1.keys(), mydic_2.keys())
	for ID in intersection:
		out = mydic_1[ID] + mydic_2[ID]
		outdic[ID]=out
	return 	outdic		
	
	
def build_list(IDs, c, file):
	"""
	returns a subset of IDs that are in column c of file.
	returns a list of (ID,annotation) 
	"""
	mylist = []
	f = open(file,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split('\t');
			if c  < len(sline):
				if sline[c] in IDs:
					item = sline[:c] + sline [c+1:]
					mylist.append((sline[c], item))
	f.close()
	return mylist

def remove_redundancy(mylist):
	"""
	list item: (key, annotation)
	If multiple elements have the same key, only the first element is retained.
	"""
	unique_list=[];
	if mylist !=[]:
		if Utility_extended.is_listT_sorted(mylist) != 1:
			mylist = sorted(mylist, key=itemgetter(0));
		unique_list.append(mylist[0]);
		for i in range(1, len(mylist)):
			if mylist[i][0] != mylist[i-1][0]:
				unique_list.append(mylist[i]);
	return unique_list;
			
def overlap (list1, list2):
	mylist = []
	for item in list1:
		if item in list2:
			mylist.append(item)
	return mylist