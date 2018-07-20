#!/usr/bin/env python

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")
	
class Records(object):
	"""
	{id:{attribute_name: attribute_value}}
	"""
	def __init__(self, records):
		self.records = copy.deepcopy(records)
		self.ids = self.records.keys()
		self.num_records = len(self.records.keys())
		self.num_attributes = len((self.records[0]).keys())
		self.attributes = sorted (self.records[0]).keys())
		self.ids_by_chrom = self.separate_ids_by_chrom()
		
	def add_attributes(self, id_dic):
		"""
		When new attributes are added, only ids in the intersection are retained
		"""
		new_ids = id_dic.keys()
		myids = self.check_id_set(new_ids)
		new_records = {}
		for myid in myids:
			new_records[myid] = self.records[myid]
			new_records[myid].update(id_dic[myid])
		self.records = copy.deepcopy(new_records)
		self.ids = self.records.keys()
		self.num_records = len(self.records.keys())
		self.num_attributes = len((self.records[0]).keys())
		self.attributes = sorted (self.records[0]).keys())
		self.ids_by_chrom = self.separate_ids_by_chrom()
		del new_records
	
	def extract(self, keys, d):
		return dict((k, d[k]) for k in keys if k in d)
	
	def output(self, outfile, ids=None):
		"""
		Output entrez genes represented by [ids] in file
		"""
		myids = self.check_id_set(ids)
		myids = sorted(myids)
		
		#header line, ordered according to alphabetical 
		outf = open(outfile, 'w')
		outline = "#" + "\t".join([str(attribute) for attribute in self.attributes]) + "\n"
		outf.write(outline)
		
		#"body", ordered according to alphabetical 
		for myid in myids:
			record = self.records[myid]
			outline = "\t".join([str(record[attribute]) for attribute in self.attributes]) + "\n"
			outf.write(outline)
		outf.close()
	
	def output_pickle(self, outfilename, ids=None):
		"""
		pickle entrez genes represented by [ids] in file
		"""
		myids = self.check_id_set(ids)
		mysubset = self.subset(myids)
		output = open(outfilename, 'wb')
		pickle.dump(mysubset, output)
		print "The number of records output to %s is %d " % (outfilename, len(myids))
		output.close()	
		
	def subset(self, ids):
		# ids is a list of entrez_ids
		# return a records dic
		myset=self.extract(ids, self.records)
		return myset
	
	def separate_ids_by_chrom(self, ids=None):
		"""
		Generate a dic {chrom:[entrez_ids]}
		"""
		myids = self.check_id_set(ids)
		ids_by_chrom = {}
		for myid in myids:
			record = self.records[myid]
			if record["chrom"] in self.chroms:
				if record["chrom"] not in ids_by_chrom.keys():
					ids_by_chrom[record["chrom"]] = []
				ids_by_chrom[record["chrom"]].append(myid)
		return ids_by_chrom
	
	def ids_by_chrom(self, chrom, ids=None):
		# myids is a dic of ids keyed by chrom
		# return a records dic
		myids = (self.separate_ids_by_chrom(ids))[chrom]
		return self.subset(myids)
	
	def get_strand_specific_ids(self, strand, ids=None):
		"""
		returns [ids]
		"""
		myids = self.check_id_set(ids)
		mygenes=[]
		assert (strand == "+" or strand =="-")
		for myid in myids:
			record = self.records[myid]
			if record["strand"] == strand:
				mygenes.append(myid)
		return mygenes
