#!/usr/bin/env python
########################################################################
## 04/01/2022
## By Xiang Li,
## lux@gwu.edu
########################################################################
#Usage 
#python findOverlaps.py -a <Bed_File1.bed> -b <Bed_File2.bed> -o <output>
########################################################################

import sys, csv, multiprocessing
from optparse import OptionParser
from itertools import groupby

########################################################################
## FUNCTIONS
### FUNCTION

def Bed_Parser(_Path):
	""" Read from input path and return as a dictionary with format: 
	[{'idx': 0, 'chr': 'chr1', 'start': '0', 'end': '5'}
	...
	"""
	df = []
	with open(_Path) as file:
		idx = 0
		for line in file:
			(chrom, start, end) = line.split()[0:3]
			df.append(({'idx': idx, 'chr': chrom, 'start': start, 'end': end}))
			idx += 1
	df = sorted(df, key = key_func)
	
	df_groups = {}
	df_keys = []
	for key, value in groupby(df, key_func):
		df_groups[key] = list(value)
		df_keys.append(key)
	return df_keys, df_groups


def key_func(k):
	## fuction for sorted by chr 
	return k['chr']


def key_func_start(k):
	## key sorted by start
	return int(k['start'])


def Bed_Intersect(_df1, _df2):
	"""Find intersection with time complexity NLogN
	Input two list with following format
	[{'idx': 0, 'chr': 'chr1', 'start': '0', 'end': '5'},
	 ...
	 {'idx': 4, 'chr': 'chr1', 'start': '30', 'end': '40'}]
	"""
	idx1, idx2 = 0, 0
	df_Intersection = []

	while idx1 < len(_df1):
		while idx2 < len(_df2):
			if (idx1 >= len(_df1)):
				break
			reg1, reg2 = _df1[idx1], _df2[idx2]
			## intersection condition
			S1, E1 = int(reg1['start']), int(reg1['end'])
			S2, E2 = int(reg2['start']), int(reg2['end'])
			if E2 < S1:
				idx2 += 1
			elif S2 > E1:
				idx1 += 1
			else:
				## Save intersection for overlap case
				df_Intersection.append([reg1['chr'], max(S1,S2), min(E1, E2)])
				idx1 += 1
		idx1 += 1
	return df_Intersection


def Write_File(_df_output, _file_name):
	with open(_file_name, 'w', newline='') as file:

		spamwriter = csv.writer(file, delimiter='\t',
								quotechar='|', quoting=csv.QUOTE_MINIMAL)
		spamwriter.writerow(["#chr", "start", "end"])
		for x in _df_output:
			spamwriter.writerow(x)
	return None

### End of FUNCTIONs
########################################################################


def main(argv):
	desc="This tool akes 2 BED files as input and reports regions of overlaps in between the 2 files. Input Format should be: #chr	bin1	bin2	featureA"
	parser = OptionParser(description=desc)
	parser.add_option("-a", "--input1", action="store", type="string",
			dest="input_file1", help="Path to first input file in Bed format", metavar="<file>")
	parser.add_option("-b", "--input2", action="store", type="string",
			dest="input_file2", help="Path to second input file in Bed format", metavar="<file>")
	parser.add_option("-o", "--output_name", action="store", type="string",
			dest="output_name", help="Name of File for Intersection Output.", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 3:
		parser.print_help()
		sys.exit(1)

## Parse Bed file format 
	df1_keys, df1_groups = Bed_Parser(opt.input_file1)
	df2_keys, df2_groups = Bed_Parser(opt.input_file2)

## only do intersection for shared chrmosome
## in case of large input datasets, can go parallel by chr
	df_shared_chr = list(set(df1_keys) & set(df2_keys))
	
	num_cpu = multiprocessing.cpu_count()
	N_thread = min(len(df_shared_chr), num_cpu)
	
	## If threads of server is not enough to run chrs at the same time. 
	## Choose CPU core as max.
	pool = multiprocessing.Pool(N_thread)
	Output_intersect = []
	
	for chrom in df_shared_chr:
		## Sorted bed file in each chromosome.
		df1 = sorted(df1_groups[chrom], key=key_func_start)
		df2 = sorted(df2_groups[chrom], key=key_func_start)
		result = pool.apply_async(Bed_Intersect, args=(df1, df2)).get()
		Output_intersect.extend(result)
	pool.close()
	pool.join()
## Sort Output File by Chromsome
	df_output = sorted(Output_intersect, key=lambda chrom: chrom[0])

	Write_File(df_output, opt.output_name+'.bed')
	print ("Output File can be found: " + opt.output_name+'.bed')

if __name__ == "__main__":
	main(sys.argv)
