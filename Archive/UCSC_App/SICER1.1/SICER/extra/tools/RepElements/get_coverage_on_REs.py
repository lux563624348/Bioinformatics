#!/usr/bin/env python

"""
Calculate the coverage of histone modifications on repetitive elements. The motivation is to be able to distinguish peakish histone mod from domainish histone mod. In particular, H3K9me3 on L1 is peakish, and on RLTR4 is spread out. 

Input:  
bedfile: island bed file
RE-tree is built by Characterize_RepElements.py
If we are only interested in a subset of REs, one should generate an appropriate RE tree.

Currently the data structure is {reClass:{reFamily:{reName:{id:value}}}}.


"""


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time

try:
	import cPickle as pickle
except:
	import pickle

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')


sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements')

import UCSC_revised
import GenomeData
import BED_revised
import Utility_extended
import RepElements
import get_read_count_on_REs

plus = re.compile("\+");
minus = re.compile("\-");


def get_coverage(re_file_name,  chrom, chrom_length, islands, upstream, downstream, min_re_length = 10):
	"""
	Find the coverage for each re instances  
	returns {id:value}
	"""
	known_repelements = RepElements.KnownRepElements.initiate_from_file([chrom], re_file_name)
	#Get the regions defined by the REs
	regions = [] #(start, end, myelement.id)
	for myid in known_repelements.rep_elements.keys():
		myelement = known_repelements.rep_elements[myid]
		#No matter whether positive or negative, genoStart < genoEnd
		if plus.match(myelement.strand):
			start = max(myelement.genoStart - upstream, 0)
			end = min(myelement.genoEnd + upstream, chrom_length)
		elif minus.match(myelement.strand):
			start = max(myelement.genoStart - downstream, 0)
			end = min(myelement.genoEnd + upstream, chrom_length)
		else:
			print myelement
			print "strand not recognized"
			exit(1)
		regions.append((start, end, myelement.id))

	#{region:coverage}, region:(start, end, myelement.id)
	regions_w_coverage = Utility_extended.find_coverage_by_islands_on_regions(regions, islands)
	#change to {id:coverage}
	coverage_dic = {}
	for myregion in regions_w_coverage.keys():
		start = myregion[0]
		end = myregion[1]
		length = end - start + 1
		if length >= min_re_length:
			myid = myregion[2]
			coverage_dic[myid] = regions_w_coverage[myregion]

	return coverage_dic


def test(summary, repClass, repFamily, repName, outputfile):
	"""
	summary: {repClass:{repFamily:{repName:{repID:value}}}}
	output the coverage data of particular (repClass, repFamily, repName)
	"""
	of = open (outputfile, 'w')
	for myid in summary[repClass][repFamily][repName].keys():
		coverage  = summary[repClass][repFamily][repName][myid]
		oline = str(myid) + "\t" + str(coverage) + "\n"
		of.write(oline)
	of.close()

def breakdown_and_output(coverage, name):
	"""
	breakdown coverage to the reName level and store the coverage info for each name separately
	
	Input:
	coverage: {reClass:{reFamily:{reName:{id:value}}}}
	
	Output:
	files storing {id:value}
	"""
	for reClass in coverage.keys():
		for reFamily in coverage[reClass].keys():
			for reName in coverage[reClass][reFamily].keys():
				file_name = name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				coverage_dic = coverage[reClass][reFamily][reName]
				output = open(file_name, 'wb')
				pickle.dump(coverage_dic, output)
				output.close()					

def main(argv):
	parser = OptionParser()

	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="Island bed file")
	parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension", help="upstream extension from start", metavar="<int>")
	parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension", help="downstream extension from end",  metavar="<int>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-n", "--feature_name", action="store", type="string", dest="feature_name",help="name of the library", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
		parser.print_help()
		sys.exit(1)

	startTime = time.time()
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);


	#load the RE tree to get the RE file names
	current_dir = os.getcwd()
	os.chdir("/home/data/mm9/Annotation/RepElements")
	re_tree = pickle.load(open("mm9_re_tree.pkl", 'rb'))
	(numb_classes, numb_families, numb_names) = get_read_count_on_REs.numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	os.chdir(current_dir)

	#Prepare the summary
	coverage = {}
	min_re_length = 10
	for reClass in re_tree.keys():
		coverage [reClass] = {}
		for reFamily in re_tree[reClass].keys():
			coverage[reClass][reFamily] = {}
			for reName in re_tree[reClass][reFamily]:
				coverage[reClass][reFamily][reName] = {}

	#Load in the island BED_GRAPH file: chrom, start, end, value
	#island.bed_vals: {chrom:[BED_GRAPH class instances]}
	islands = BED_revised.BED(opt.species, opt.bedfile, "BED_GRAPH");

	# chdir to load in the individual RE annotation file
	current_dir = os.getcwd()
	os.chdir("/home/data/mm9/Annotation/RepElements/re_by_name")			
	#cycle through chrom
	for chrom in chroms:
		chrom_length = chrom_lengths[chrom]
		print chrom
		# [(start, end, value)]
		chrom_specific_islands =[ (item.start, item.end, item.value) for item in islands.bed_vals[chrom]]
		i = 0
		for reClass in re_tree.keys():
			for reFamily in re_tree[reClass].keys():
				for reName in re_tree[reClass][reFamily]:
					re_file_name = "_".join([reClass, reFamily, reName]) + ".txt"
					#{id:value}
					i += 1
					#print reClass, reFamily, reName, i
					coverage_dic = get_coverage(re_file_name, chrom, chrom_length, chrom_specific_islands, opt.upstream_extension, opt.downstream_extension, min_re_length) 
					# id is unique and updated only once, so this should be ok
					coverage[reClass][reFamily][reName].update(coverage_dic)
	os.chdir(current_dir)
	#{reClass:{reFamily:{reName:{id: value}}}
	#output_file_name = feature_name + "_on_" + "mm9_rmsk.pkl"
	#output = open(output_file_name, 'wb')
	#pickle.dump(coverage, output)
	#output.close()	

	#instead of outputing a huge one, let's output many small pieces
	breakdown_and_output(coverage, opt.feature_name)

	repClass ='LTR'
	repFamily = 'ERV1'
	repName = 'RLTR4_Mm'
	outfile_name =opt.feature_name  + "_on_" +  "_".join([repClass, repFamily, repName]) + ".dat"
	test(coverage, repClass, repFamily, repName, outfile_name)
	print repClass, repFamily, repName

	repClass ='LINE'
	repFamily = 'L1'
	repName = 'L1Md_Gf'
	outfile_name =opt.feature_name  + "_on_" +  "_".join([repClass, repFamily, repName]) + ".dat"
	test(coverage, repClass, repFamily, repName, outfile_name)
	print repClass, repFamily, repName


	print "it took", time.time() - startTime, "seconds."


if __name__ == "__main__":
	main(sys.argv)
