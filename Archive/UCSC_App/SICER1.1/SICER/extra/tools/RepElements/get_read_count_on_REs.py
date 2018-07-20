#!/usr/bin/env python

"""

RE-tree is built by Characterize_RepElements.py

Calculate read count on repetitive elements
Currently the data structure is {reClass:{reFamily:{reName:{id:{feature_name:value}}}}}, the memory usage is high and the speed is quite slow. With two features, the hard drive storage space is three times the storage used by {reClass:{reFamily:{reName:{id:value}}}}. The expected should be two times. The memory problem gets even worse when features are combined, 8G is not enough for combining index 6 and index 8. We need to do something. 
Solution: store read count for each type RE separately.
This is already implemented

If we are only interested in a subset of REs, one should generate an appropriate RE tree.

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
import SeparateByChrom
import Utility_extended
import RepElements
import associate_tags_with_regions
import get_total_tag_counts

plus = re.compile("\+");
minus = re.compile("\-");

def numbers(re_tree):
	numb_classes = len(re_tree.keys())
	numb_families = sum([len(item.keys()) for item in re_tree.values()])
	numb_names = 0
	for reClass in re_tree.keys():
		for reFamily in re_tree[reClass].keys():
			numb_names += len(re_tree[reClass][reFamily])
	return (numb_classes, numb_families, numb_names)


def get_read_count(re_file_dir, re_file_name, feature_name,  chrom, chrom_length, tag_position_list, total_count, upstream, downstream, min_re_length):
	"""
	returns {id:{feature_name: value}}
	feature_name include: feature_name + "_rc", feature_name + "_rpkm"
	
	"""
	currentdir = os.getcwd()
	os.chdir(re_file_dir)
	known_repelements = RepElements.KnownRepElements.initiate_from_file([chrom], re_file_name)
	
	regions = []
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
		
	tag_list = [(element,0) for element in tag_position_list]
	read_counts = Utility_extended.get_read_counts_on_regions(tag_list, regions) #returns a list, [region, read_count]
	rc_dic = {}
	for item in read_counts:
		region = item[0]
		start = region[0]
		end = region[1]
		myid = region[2]
		if start == end:
			print chrom, myid, start, end
		if (end-start) >= min_re_length: #only include those with length >= min_re_length
			rc = item[1]
			rpkm = rc / ((total_count)/1000000.0)
			rpkm = rpkm / ((end-start)/1000.0)
			rc_dic[myid] = {}
			rc_dic[myid][feature_name + "_rc"] = rc
			rc_dic[myid][feature_name + "_rpkm"] = rpkm
	os.chdir(currentdir)
	return rc_dic

def characteristics(feature):
	"""
	data structure {reClass:{reFamily:{reName:{id:{feature_name:value}}}}}
	"""
	numb_classes = len(feature.keys())
	numb_families = sum([len(item.keys()) for item in feature.values()])
	numb_names = 0
	numb_ids = 0

	for reClass in feature.keys():
		for reFamily in feature[reClass].keys():
			numb_names += len(feature[reClass][reFamily].keys())
			for reName in feature[reClass][reFamily].keys():
				numb_ids += len(feature[reClass][reFamily][reName].keys())
	
	my_reClass = feature.keys()[0]
	my_reFamily = feature[my_reClass].keys()[0]
	my_reName = feature[my_reClass][my_reFamily].keys()[0]
	my_id = feature[my_reClass][my_reFamily][my_reName].keys()[0]
	features = feature[my_reClass][my_reFamily][my_reName][my_id].keys()
	
	return (numb_classes, numb_families, numb_names, numb_ids, features)
		

def test(summary, repClass, repFamily, repName, outputfile):
	"""
	summary: {repClass:{repFamily:{repName:{repID:{feature:value}}}}}
	output the rc data of particular (repClass, repFamily, repName)
	"""
	of = open (outputfile, 'w')
	myres = summary[repClass][repFamily][repName]
	myids = myres.keys()
	myid = myids[0]
	myfeatures = myres[myid].keys()
	myfeatures = sorted (myfeatures)
	oline = ("\t").join(myfeatures) + "\n"
	of.write(oline)
	for myid in summary[repClass][repFamily][repName].keys():
		re  = summary[repClass][repFamily][repName][myid]
		oline = str(myid)	
		for feature in myfeatures:
			oline +=  "\t" + str(re[feature])
		oline += "\n"
		of.write(oline)
	of.close()

def breakdown_and_output(rc_dic, name):
	"""
	rc_dic: {reClass:{reFamily:{reName:{id:{feature_name:value}}}}}
	breakdown rc_dic to the reName level and store the rc info for each name separately
	After breaking down: {id:{feature_name:value}}
	"""
	for reClass in rc_dic.keys():
		for reFamily in rc_dic[reClass].keys():
			for reName in rc_dic[reClass][reFamily].keys():
				file_name = name + "_on_" + "_".join([reClass, reFamily, reName]) + ".pkl"
				res = rc_dic[reClass][reFamily][reName]
				output = open(file_name, 'wb')
				pickle.dump(res, output)
				output.close()					
	
	
def main(argv):
	parser = OptionParser()
	
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-t", "--RE_tree_pickle_file", action="store", type="string", dest="RE_Tree", metavar="<file>", help="file with RE tree in pickle format")
	parser.add_option("-l", "--RE_annotation_file_location", action="store", type="string", dest="RE_file_location", metavar="<file>", help="location of RE files named in repClass_repFamily_repName.txt")
	parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension", help="upstream extension from start", metavar="<int>")
	parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension", help="downstream extension from end",  metavar="<int>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-n", "--feature_name", action="store", type="string", dest="feature_name",help="name of the library", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	total_count = get_total_tag_counts.get_total_tag_counts(opt.bedfile);
	
	#Separate_by_chrom on bedfile
	lib_name = (opt.bedfile).split('/')[-1] # remove directory
	suffix = lib_name.split('.')[-1] # txt
	lib_name = lib_name.split('.')[0] 
	extension = "-" + lib_name +'.' + suffix +"1"
	if Utility_extended.fileExists(opt.bedfile):
		if Utility_extended.chrom_files_exist(chroms, extension) != 1:
			SeparateByChrom.separateByChrom(chroms, opt.bedfile,  extension)
	else:
		print bedfile, " is not found";
		sys.exit(1)
	
	
	
	#load the RE tree to get the RE file names
	re_tree = pickle.load(open(opt.RE_Tree, 'rb'))
	(numb_classes, numb_families, numb_names) = numbers(re_tree)
	print "There are %d classes, %d family, and %d names." %(numb_classes, numb_families, numb_names)
	
	
	#Prepare the summary
	read_counts = {}
	for reClass in re_tree.keys():
		read_counts [reClass] = {}
		for reFamily in re_tree[reClass].keys():
			read_counts[reClass][reFamily] = {}
			for reName in re_tree[reClass][reFamily]:
				read_counts[reClass][reFamily][reName] = {}
				
	#cycle through chrom
	for chrom in chroms:
		print chrom
		chrom_length = chrom_lengths[chrom]
		chrombed = chrom + extension
		if Utility_extended.fileExists(chrombed):		
			# load in each read and shift
			tag_position_list = []
			inf = open(chrombed,'r')
			for line in inf:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_position_list.append(associate_tags_with_regions.tag_position(sline, opt.fragment_size))
			inf.close()	
			if not Utility_extended.is_list_sorted(tag_position_list):
				tag_position_list.sort() #[tag_positions]
	
		min_re_length = 10
		for reClass in re_tree.keys():
			for reFamily in re_tree[reClass].keys():
				for reName in re_tree[reClass][reFamily]:
					re_file_name = "_".join([reClass, reFamily, reName]) + ".txt"
					#{id:{feature_name:value}}
					rc_dic = get_read_count(opt.RE_file_location, re_file_name, opt.feature_name, chrom, chrom_length, tag_position_list, total_count,  opt.upstream_extension, opt.downstream_extension, min_re_length) 
					# id is unique and updated only once, so this should be ok
					read_counts[reClass][reFamily][reName].update(rc_dic)
	
	#{reClass:{reFamily:{reName:{id:feature_name, value}}}}
	#feature_name include: feature_name + "_rc", feature_name + "_rpkm"
	#output_file_name = feature_name + "_on_" + "mm9_rmsk.pkl"
	#output = open(output_file_name, 'wb')
	#pickle.dump(read_counts, output)
	#output.close()	
	
	#instead of outputing a huge one, let's output many small pieces
	breakdown_and_output(read_counts, opt.feature_name)
	
	repClass ='LTR'
	repFamily = 'ERV1'
	repName = 'RLTR4_Mm'
	outfile_name = lib_name + "_on_" +  "_".join([repClass, repFamily, repName]) + ".dat"
	test(read_counts, repClass, repFamily, repName, outfile_name)
	
	SeparateByChrom.cleanup(chroms, extension)
	
	print "it took", time.time() - startTime, "seconds."

	
if __name__ == "__main__":
	main(sys.argv)
