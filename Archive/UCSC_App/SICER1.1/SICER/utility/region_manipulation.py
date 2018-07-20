import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import Utility
from GenomeData import *
import SeparateByChrom

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def are_regions_sorted(region_list):
	"""
	Check if sorted in ascending order.
	input is a list of BED with chrom, start, end and value.
	output: sorted =1 or 0
	"""
	sorted = 1;
	for index in range(0, len(region_list)-1):
		if region_list[index].start> region_list[index+1].start:
			sorted = 0;
	return sorted;

def are_regions_overlapping(region_list):
	"""
	Check if regions are overlapping
	"""
	if are_regions_sorted(region_list) == 0:
		region_list.sort(key = operator.attrgetter('start'))
	
	for i in xrange(len(region_list)-1):
		if region_list[i].end >= region_list[i+1].start: # assume closed end
			return 1
	return 0

def union_islands(islandlist):
	"""
	The islandlist MUST be pre-sorted according to start!!!
	"""
	start_list =[]
	end_list = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			start_list.append(current.start)
			end_list.append(current.end)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	start_list.append(current.start)
	end_list.append(current.end)
	return start_list, end_list