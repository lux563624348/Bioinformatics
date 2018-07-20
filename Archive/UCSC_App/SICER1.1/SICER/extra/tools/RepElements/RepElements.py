#!/usr/bin/env python
# Author :  Weiqun Peng


"""
This module contains the classes to deal with UCSC repetitive
element files.
"""

import re, os, shutil, time, sys, operator
from math import *
from string import *
from optparse import OptionParser
import copy
from operator import itemgetter
try:
   import cPickle as pickle
except:
   import pickle

#import numpy
#import matplotlib.pyplot as plt
#import matplotlib

sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import GenomeData 
import UCSC_revised
import Utility_extended

plus = re.compile('\+');
minus = re.compile('\-');

RepElementsError = "Error in RepElement class";



class RepElement:
	"""
	Class for storing repetitive element information.
	
	The "repStart and repEnd" coordinates represent the region of the query sequence (RepBase sequence aka repName) involved in the alignment to genomic. 
	
	Of the three coordinates involved (repStart, repEnd, repLeft), repEnd can be thought of as the "middle" coordinate of the three.  It is always a positive integer regardless of the strand the repeat element hits.  It represents the end coordinate of the matching part of the repeat element in "repeat element coordinates".   The "beginning" coordinate of the matching part of a repeat element is repStart (for + oriented hits) or repLeft (for - oriented hits).

	For + oriented hits, repStart & repEnd are the proper coordinates of the matching part of the element.  RepLeft in this case is a numerical value which may be used (via the equation repEnd-repLeft) to obtain the size of the repeat element ("Left" in the sense of "repeat remaining" unaligned).

	For - oriented hits, repLeft & repEnd are the proper coordinates of the matching part of the element.  RepLeft will be upstream from (smaller than) repEnd for these neg-strand alignments.  RepStart in this case is a numerical value which may be used (via the equation repEnd-repStart) to obtain the size of the repeat element.
	
	coord.id = coord.repName + chrom  + str(coord.strand) + str(coord.genoStart)
	
	# calculate the fraction of match in concensus rep sequence
	if re.match(plus, coord.strand):
		coord.quality = (coord.repEnd - coord.repStart)/(float(coord.repEnd - coord.repLeft));
	elif re.match(minus, coord.strand):
		coord.quality = (coord.repEnd - coord.repLeft)/(float(coord.repEnd -coord.repStart));
	
	#Calculate length of match in genome
	coord.length = coord.genoEnd - coord.genoStart;
	"""
	def __init__(self, bin, swScore, milliDiv, milliDel, milliIns, chrom,
				genoStart, genoEnd, genoLeft, strand,
				repName, repClass, repFamily, repStart, repEnd,
				repLeft, id, age, quality, length):

		self.bin = bin; #0
		self.swScore = swScore;
		self.milliDiv = milliDiv;
		self.milliDel = milliDel;
		self.milliIns = milliIns;
		self.chrom = chrom #5
		self.genoStart = genoStart;
		self.genoEnd = genoEnd;
		self.genoLeft = genoLeft;
		self.strand = strand; #9
		self.repName = repName;
		self.repClass = repClass;
		self.repFamily = repFamily;
		self.repStart = repStart; #13
		self.repEnd = repEnd;#14
		self.repLeft = repLeft;
		self.id = id; #16
		self.age = age;
		self.quality = quality;
		self.length = length;

	def copy(self, a):
		self.bin = a.bin;
		self.swScore = a.swScore;
		self.milliDiv = a.milliDiv;
		self.milliDel = a.milliDel;
		self.milliIns = a.milliIns;
		self.chrom = a.chrom
		self.genoStart = a.genoStart;
		self.genoEnd = a.genoEnd;
		self.genoLeft = a.genoLeft;
		self.strand = a.strand;
		self.repName = a.repName;
		self.repClass = a.repClass;
		self.repFamily = a.repFamily;
		self.repStart = a.repStart;
		self.repEnd = a.repEnd;
		self.repLeft = a.repLeft;
		self.id = a.id;
		self.age = a.age;
		self.quality = a.quality;
		self.length = a.length;

	def __setitem__(self, bin, swScore, milliDiv, milliDel, milliIns, chrom,
				genoStart, genoEnd, genoLeft, strand,
				repName, repClass, repFamily, repStart, repEnd,
				repLeft, id, age, quality, length):

		self.bin = bin;
		self.swScore = swScore;
		self.milliDiv = milliDiv;
		self.milliDel = milliDel;
		self.milliIns = milliIns;
		self.chrom = chrom
		self.genoStart = genoStart;
		self.genoEnd = genoEnd;
		self.genoLeft = genoLeft;
		self.strand = strand;
		self.repName = repName;
		self.repClass = repClass;
		self.repFamily = repFamily;
		self.repStart = repStart;
		self.repEnd = repEnd;
		self.repLeft = repLeft;
		self.id = id;
		self.age = age;
		self.quality = quality;
		self.length = length;
		
	def __str__(self):
		"""
		Prints out all important given information.  Information used for
		sorting or filtering not given (quality, age, length of sequence)
		"""
		outstring = str(self.bin) + " " + str(self.swScore) + " " + \
					str(self.milliDiv) + " " + str(self.milliDel) + " " + \
					str(self.milliIns) + " " + self.chrom + " "  + \
					str(self.genoStart) + " " + str(self.genoEnd) + " " + \
					str(self.genoLeft) + " " + self.strand + " " + \
					self.repName + " " + self.repClass + " " + \
					self.repFamily + " " + str(self.repStart) + " " + \
					str(self.repEnd) + " " + str(self.repLeft) + " " + \
					str(self.id) + '\n';
		try:
			return outstring;
		except:
			sys.stderr.write("No UCSC repetitive element information for %s\n" % self)
			return ''

#--------------------------
#--------------------------

class KnownRepElements:
	"""
	Class to read in repetitive element files and organize information.
	{id:RepElement}
	"""

	def __init__(self, dic={}):
		"""
		Two ways of initiation:
			Reads in a repElement file and builds a dictionary with chromosomes as keys and all other information appended in a list as values.
			copy of a made dictionary. This is very helpful for subset operation
			
		
		"""
		self.rep_elements = copy.deepcopy(dic)
		self.chroms = self.get_chroms()
		
		self.repClass_statistics = self.get_repClass_statistics()
		self.repClass_names = self.repClass_statistics.keys()
		
		self.repFamily_statistics = self.get_repFamily_statistics()
		self.repFamily_names = self.repFamily_statistics.keys()
		
		self.repName_statistics = self.get_repName_statistics()
		self.repName_names = self.repName_statistics.keys()
		
		self.number = len(self.rep_elements.keys())
		
		#print "There are %d classes, %d families, and %d rep elements" % (len(self.repClass_names), len(self.repFamily_names), len(self.repName_names))
		#print self.repClass_names 
		#print self.repFamily_names
		#print self.repName_names
		#mylist = self.print_statistics(self.repName_statistics)
		#for item in mylist:
			#print item[0], item[1]
		
	@classmethod
	def initiate_from_file(cls, chroms, file): # use chroms to remove contribution from exotic chroms
		"""
		Initiate from file 
		chroms are used to remove things on exotic chroms 
		
		
		The "repStart and repEnd" coordinates represent the region of the query sequence (RepBase sequence aka repName) involved in the alignment to genomic. 
	
		Of the three coordinates involved (repStart, repEnd, repLeft), repEnd can be thought of as the "middle" coordinate of the three.  It is always a positive integer regardless of the strand the repeat element hits.  It represents the end coordinate of the matching part of the repeat element in "repeat element coordinates".   The "beginning" coordinate of the matching part of a repeat element is repStart (for + oriented hits) or repLeft (for - oriented hits).

		For + oriented hits, repStart & repEnd are the proper coordinates of the matching part of the element.  RepLeft in this case is a numerical value which may be used (via the equation repEnd-repLeft) to obtain the size of the repeat element ("Left" in the sense of "repeat remaining" unaligned).

		For - oriented hits, repLeft & repEnd are the proper coordinates of the matching part of the element.  RepLeft will be upstream from (smaller than) repEnd for these neg-strand alignments.  RepStart in this case is a numerical value which may be used (via the equation repEnd-repStart) to obtain the size of the repeat element.
	
		age is defined as (atoi(sline[2]) + atoi(sline[3]) + atoi(sline[4])): milliDiv, milliDel, milliIns
		
		"""
		rep_elements = {}
		infile = open(file);
		for line in infile:
			if not re.match("#", line):
				line = line.strip();   #"""Strip white space off line""" 
				sline = line.split();  #"""Split line into individual strings (fields)"""
				""" Check to make sure this chromosome is declared in the dictionary."""
				chrom = sline[5]
				if chrom in chroms: # Remove exotic chroms
					coord = RepElement( atoi(sline[0]), atoi(sline[1]), atoi(sline[2]), atoi(sline[3]), atoi(sline[4]), sline[5], atoi(sline[6]), atoi(sline[7]), sline[8], sline[9], sline[10], sline[11], sline[12],
										atoi(sline[13]), atoi( sline[14]), atoi(sline[15]), sline[16],
										(atoi(sline[2]) + atoi(sline[3]) + atoi(sline[4])), #age
										0, 0);
					#The original id is not really an id, we here redefine it. 
					coord.id = coord.repName + chrom  + str(coord.strand)  + str(coord.genoStart)
					#Check to see if the strand is positive or negative, and calculate quality accordingly.
					if re.match(plus, coord.strand):
						coord.quality = (coord.repEnd - coord.repStart)/(float(coord.repEnd - coord.repLeft));
					elif re.match(minus, coord.strand):
						coord.quality = (coord.repEnd - coord.repLeft)/(float(coord.repEnd -coord.repStart));
					"""Calculate length of match on genome"""
					coord.length = coord.genoEnd - coord.genoStart;
					"""Add coords to the dictionary under the indicated id"""
					rep_elements[coord.id] = coord
					
		return cls(rep_elements)
	
	def output(self, outfile, ids=None):
		universe = self.check_id_set (ids)
		outf = open(outfile, 'w')
		for myid in universe:
			element = self.rep_elements[myid]
			out.write(str(element))
		outf.close()
	
	def output_pickle(self, outfilename, ids=None):
		"""
		pickle re represented by [ids] in file
		"""
		universe = self.check_id_set (ids)
		mysubset = self.subset(universe)
		output = open(outfilename, 'wb')
		pickle.dump(mysubset, output)
		print "The number of elements output to %s is %d " % (outfilename, len(universe))
		output.close()	
	
	def get_rep_tree(self, ids=None):
		"""
		Build the tree structure
		"""
		self.mytree = {}
		universe = self.check_id_set (ids)
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repClass not in self.mytree.keys():
				self.mytree[element.repClass] = {}
				self.mytree[element.repClass][element.repFamily] = []
				self.mytree[element.repClass][element.repFamily].append(element.repName)
			else:
				if element.repFamily not in self.mytree[element.repClass]:
					self.mytree[element.repClass][element.repFamily] = []
					self.mytree[element.repClass][element.repFamily].append(element.repName)
				else:
					if element.repName not in self.mytree[element.repClass][element.repFamily]:
						self.mytree[element.repClass][element.repFamily].append(element.repName)
		for myclass in self.mytree.keys():
			print myclass
			for myfamily in self.mytree[myclass].keys():
				print "\t", myfamily
				for myname in self.mytree[myclass][myfamily]:
					print "\t\t", myname

	
	def get_repClass_statistics(self, ids=None):
		"""
		if None, do it for all
		Return statistics for repClasses: {repClass:count}
		"""
		universe = self.check_id_set(ids)
		repClasses = {}
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repClass not in repClasses.keys():
				repClasses[element.repClass] = 1
			else:
				repClasses[element.repClass] += 1
		return repClasses
	
	def find_class_members(self, class_name, ids=None):
		"""
		Finds all repetitive elements carrying a certain name 
		return [ids]
		"""
		
		universe = self.check_id_set (ids)
		
		myid = []; # Dict. to  store all elem. of specified repClass.		
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repClass == class_name:
				myid.append(myid)
		return myid
					
	def get_repFamily_statistics(self, ids=None):
		"""
		Return statistics for repFamilies: {repFamilyNames:count}
		"""
		universe = self.check_id_set (ids)
			
		repFamilies = {}
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repFamily not in repFamilies.keys():
				repFamilies[element.repFamily] = 1
			else:
				repFamilies[element.repFamily] += 1
		return repFamilies
	
	def find_family_members(self, family_name, ids=None):
		"""
		Finds all repetitive elements in a certain family and inserts
		returns ids
		"""
		universe = self.check_id_set (ids)

		myid = []; # Dict. to  store all elem. of specified family.
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repFamily == family_name:
				myid.append(myid)
		return myid
	
	
	def get_repName_statistics(self, ids=None):
		"""
		Return statistics for repNames: {repName:count}
		"""
		universe = self.check_id_set (ids)
			
		repNames = {}
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repName not in repNames.keys():
				repNames[element.repName] = 1
			else:
				repNames[element.repName] += 1
		return repNames
	
	def find_name_members(self, name, ids=None):
		"""
		Finds all repetitive elements carrying a certain name 
		Return [ids]
		"""
		universe = self.check_id_set (ids)
					
		myids = []; 
		for myid in universe:
			element = self.rep_elements[myid]
			if element.repName == name:
				myids.append(myid)
		return myids

	def print_statistics(self, statistics):
		"""
		print out the three statistics get_repClass_statistics, get_repFamily_statistics, get_repName_statistics
		"""
		mylist = [(myid, statistics[myid]) for myid in statistics.keys()]
		mylist = sorted(mylist, key = itemgetter(1), reverse=True)
		return mylist
	
	
	def get_age_distribution(self, ids=None):
		universe = self.check_id_set (ids)
		my_list = [self.rep_elements[myid].age for myid in universe] 
		return my_list
		
	def average_age(self, ids=None):
		"""
		Finds average age 
		"""
		my_list = self.get_age_distribution(ids) 
		avg_age = sum(my_list) / float(len(my_list));
		return avg_age;
	
	def sort_by_age(self, ids=None):
		"""
		Sorts by age. Age is a value that is the sum
		of milliDiv, milliDel, and milliIns.
		
		returns the ordered [ids] 
		"""
		universe = self.check_id_set (ids)		
		my_list = [self.rep_elements[myid] for myid in universe] 
		my_list.sort(key=operator.attrgetter('age'));
		return [element.id for element in my_list]
	
	def get_quality_distribution(self, ids=None):
		universe = self.check_id_set (ids)
		my_list = [self.rep_elements[myid].quality for myid in universe] 
		return my_list
		
	def average_quality(self, ids=None):
		"""
		Finds average quality 
		"""
		my_list = self.get_quality_distribution(ids) 
		avg_quality = sum(my_list) / float(len(my_list));
		return avg_quality;
	
	def sort_by_quality(self, ids=None):
		"""
		Sorts by the quality of the repetitive element.
		Quality is determined by a value representing the ratio of
		the match sequence to the whole repetitive element. Sorts
		highest quality to lowest.
		"""
		universe = self.check_id_set (ids)
		my_list = [self.rep_elements[myid] for myid in universe] 
		my_list.sort(key=operator.attrgetter('quality'));
		return [element.id for element in my_list]
	
	def filter_by_quality(self, quality_threshold, ids=None):
		universe = self.check_id_set (ids)
		filtered = []	
		for myid in universe:
			element = self.rep_elements[myid]
			if element.quality >= quality_threshold:
				filtered.append(myid)
		return filtered
		
	
	def filter_by_length(self, min_length=250, ids=None):
		"""
		Removes all repetitive element matches which are under
		250 base pairs.  Returns a filtered ids
		"""
		universe = self.check_id_set (ids)
		filtered = []	
		for myid in universe:
			element = self.rep_elements[myid]
			if int(element.length) >= min_length:
				filtered.append(myid);
		return filtered
		
	def get_length_distribution(self, ids=None):
		universe = self.check_id_set (ids)
		my_list = [self.rep_elements[myid].length for myid in universe] 
		return my_list

	def average_length(self, ids=None):
		"""
		Finds average length 
		"""
		my_list = self.get_length_distribution(ids)
		avg_list = sum(my_list) / float(len(my_list));
		return avg_list
	
	def get_rep_length_distribution(self, ids=None):
		my_list = []
		universe = self.check_id_set (ids)
		for myid in universe:
			element = self.rep_elements[myid]
			if re.match(plus, element.strand):
						rep_length = element.repEnd - element.repLeft
			elif re.match(minus, element.strand):
						rep_length = element.repEnd - element.repStart
		my_list.append(rep_length)
		return my_list
	
	
	def get_chroms(self, ids=None):
		universe = self.check_id_set (ids)
		chroms = list(set([self.rep_elements[myid].chrom for myid in universe]))
		return chroms
		
	def subset(self, ids):
		# ids is a list of repelement_ids
		universe = self.check_id_set (ids)
		# return a rep_elements dic
		mydic = dict((k, self.rep_elements[k]) for k in universe if k in self.rep_elements)
		return mydic
		
	
	def separate_by_chrom(self, ids=None):
		"""
		Generate a dic {chrom:[myids]}
		"""
		universe = self.check_id_set (ids)
		
		repe_ids_by_chrom = {}
		for myid in universe:
			element = self.rep_elements[myid]
			if element.chrom not in repe_ids_by_chrom.keys():
				repe_ids_by_chrom[element.chrom] = []
			repe_ids_by_chrom[element.chrom].append(myid)
		return repe_ids_by_chrom
	
	def subset_by_chrom(self, chrom):
		# myid is a dic of myids keyed by chrom
		# return a entrez_elements dic
		if chrom in self.chroms:
			return self.subset(self.repe_ids_by_chrom[chrom])
		else:
			return {}
			
	def get_strand_specific_elements(self, ids, strand):
		"""
		returns [myids]
		"""
		myids = self.check_id_set(ids)
		
		myelements=[]
		assert (strand == "+" or strand =="-")
		for id in myids:
			element = self.rep_elements[id]
			if element.strand == strand:
				myelements.append(id)
		return myelements
	
	def check_id_set(self, ids=None):
		"""
		check whether [ids] is a subset of self.rep_elements
		if not, return intersection
		"""
		if ids is None:
			universe = self.rep_elements.keys()
		else:
			universe = list(set(ids) & set(self.rep_elements.keys()))
			if len(universe) < len(ids):
				print "%d ids are thrown out" % ( len(ids) - len(universe))
		return universe
	
	
	def __del__(self):
		"""
		Delete, delete;
		"""
		self.rep_elements.clear()

	def __contains__(self, item):
		"""
		Returns  mapping iterator
		"""
		return self.rep_elements.has_key(item)

	def __iter__(self):
		"""
		Returns mapping iterator
		"""
		return self.rep_elements.iterkeys()

	def __len__(self):
		"""
		Returns number of coords
		"""
		return len(self.rep_elements.keys())
#-----------------------------------------
#-----------------------------------------

def main(argv):

	"""Lines to parse command line arguments from the file repUCSC.sh"""
	parser = OptionParser()
	parser.add_option("-i", "--file", action="store", type="string", dest="datafile", metavar="<file>",
					help="text file in raw format")    
	parser.add_option("-o","--operation", action="store", type="string", dest="op_type", metavar="<str>",
					help="name of operation to be performed")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	rep_element_collection = KnownRepElements.initiate_from_file(chroms, opt.datafile)
	print rep_element_collection.number
	rep_element_collection.get_rep_tree()	

	#min_quality= 0.7
	
	#rep_element_collection =KnownRepElements( rep_element_collection.subset(rep_element_collection.filter_by_quality(min_quality)))
	#print rep_element_collection.number
	
	##filename = "Alu_elements_q0.7.pkl"
	##rep_element_collection.output_pickle(filename)
	
	
	#plt.clf()
	#plt.figure(1)
	
	#age_list = rep_element_collection.get_age_distribution()
	#plt.subplot(411)
	#plt.hist(age_list, bins=100, color='r', normed=True, alpha= 1)
	#plt.xlabel("age")
	#plt.ylabel("Frequency")
	
	#quality_list = rep_element_collection.get_quality_distribution()
	#plt.subplot(412)
	#plt.hist(quality_list, bins=100, color='r', normed=True)
	#plt.hist(quality_list, bins=50, color='r', normed=True, log=True)
	#plt.xlabel("quality")
	#plt.ylabel("Frequency")
	#plt.legend(loc = 'upper left')
	
	#genome_length_list = rep_element_collection.get_length_distribution()
	#plt.subplot(413)
	#plt.hist(genome_length_list, bins=100, color='r', normed=True)
	##plt.hist(quality_list, bins=50, color='r', normed=True, log=True)
	#plt.xlabel("length of match on genome")
	#plt.ylabel("Frequency")
	
	#rep_length_list = rep_element_collection.get_rep_length_distribution()
	#plt.subplot(414)
	#plt.hist(rep_length_list, bins=100, color='r', normed=True)
	##plt.hist(quality_list, bins=50, color='r', normed=True, log=True)
	#plt.xlabel("rep_element length")
	#plt.ylabel("Frequency")
	
	#plt.savefig( "Alu_distributions_" + str(min_quality) +".png", format="png")
	
	#if (opt.op_type) == "KnownRepElements":
		### Add functions from KnownRepElements class
		### above to accomplish task.
		#os.chdir('/home/shane/data/sorted_result')
		#for files in os.listdir('/home/shane/data/sorted_result'):
			#print files;
			#main_dict = KnownRepElements(files);
			#main_dict.Average_Age();
			#main_dict.Average_Length();        

	#if (opt.op_type) == "Process_RawData":
		#Process_RawData();
	
if __name__ == "__main__":
	main(sys.argv)
