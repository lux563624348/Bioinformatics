# Copyright (c) 2009 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao

#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# schonesde@mail.nih.gov).


# Modified by Weiqun Peng in 4/17/2010 to 
# 1) conform to BED convention on the ordering of columns, except for BED_GRAPH format
#      "BED format provides a flexible way to define the data lines that are displayed in an annotation track. BED lines have three required fields (chrom, start, end) and nine additional optional fields (name, score, strand, ThickStart, ThickEnd, itemRgb, blockCount, blockSizes, blockStarts). The number of fields per line must be consistent throughout any single set of data in an annotation track. The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used."
# 2) remove the restriction on the number of columns 
#



"""
This module contains all the classes to deal with BED files

	BED types:
	BED2: start, strand
	BED3: chrom, start and end
	BED6: BED3 + name + score + strand ('+' or '-')
	BED_GRAPH: chrom, start, end, value
	BED_annotated: Class to deal with bed-like file, assign the needed fields to corresponding variables, and dump the rest into a field called annotation. 
  
  Additional types:
  FIXED_STEP: the memory efficient version of BED_GRAP




"""
import re, os, shutil, time, sys
from math import *
from string import *
import GenomeData;


plus = re.compile('\+');
minus = re.compile('\-');
comment = re.compile('track|#')

bedError = "Error in BED class";



#----------------------------
#----------------------------

#####   TO BE DONE   #####

## 1 - the chromosome element could be taken out of all the
## bed classes since it is already the key to the dictionary - this
## will save some memory, but this is also sometimes useful.


## 2 - BED4 should be replaced by a BEDmulti, or something like this


## 3 - BED_PAIR should be allowed to use something other than BED2
## elements

#----------------------------
#----------------------------

class BED2:
	"""
	Class for bed lines with 2 values: start and strand, this will be
	useful for things like TSS information or tags where the only
	important information is the start and strand.
	"""
	def __init__(self, start, strand):
		self.start = start;
		self.strand = strand;

	def __set__(self, start, strand):
		self.start = start;
		self.strand = strand;

	def getCoord(self):
		outstring = str(self.start) + " " + self.strand;
		try:
			return outstring;
		except:
			sys.stderr.write("No coord information for %s\n" % self)
			return ''


#----------------------------
#----------------------------

class BED3:
	"""
	Class for bed lines with 3 values: chrom, start and end
	"""
	def __init__(self, chrom, start, end):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
	def __set__(self, chrom, start, end):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
	def getCoord(self):
		outstring = self.chrom + " " + str(self.start) + " " + str(self.end);
		try:
			return outstring;
		except:
			sys.stderr.write("No coord information for %s\n" % self)
			return ''


class BED4:
	"""
	Class for bed lines with 4 values: start, end, name, strand

	** this is essentially a memory efficient version of BED6,
	not storing the chromosome -- which will usually be in the
	key of the dictionary anyways -- or the score.
	"""
	def __init__(self, start, end, name, strand):
		self.start = start;
		self.end = end;
		self.name = name;
		self.strand = strand;
	def __set__(self, start, end, name, strand):
		self.start = start;
		self.end = end;
		self.name = name;
		self.strand = strand;


#----------------------------
#----------------------------


class BED6:
	"""
	Class for bed lines with 6 values:  chrom, start, end, name, score, strand
	"""
	def __init__(self, chrom, start, end, name, score, strand):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
		self.name = name;
		self.score = score;
		self.strand = strand;
	def __set__(self, chrom, start, end, name, score, strand):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
		self.name = name;
		self.score = score;
		self.strand = strand;
	def getCoord(self):
		outstring = self.chrom + "\t" + str(self.start) + "\t" + \
					str(self.end) + "\t" + self.name + "\t" + \
					str(self.score) + "\t" + self.strand;
		try:
			return outstring;
		except:
			sys.stderr.write("No coord information for %s\n" % self)
			return ''

#----------------------------
#----------------------------

# The BED_GRAPH format does not conform to the BED format requirement of the ordering of the columns.
class BED_GRAPH:
	"""
	Class to deal with bed graph lines: chrom, start, end, value
	This emulates the wiggle format
	"""
	def __init__(self, chrom, start, end, value=0):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
		self.value = value;
	def __set__(self, chrom, start, end, value):
		self.chrom = chrom;
		self.start = start;
		self.end = end;
		self.value = value;
	def getCoord(self):
		outstring = self.chrom + " " + str(self.start) + " " + str(self.end);
		try:
			return outstring;
		except:
			sys.stderr.write("No BED coord information for %s\n" % self)
			return ''
	def getAll(self):
		outstring = self.chrom + " " + str(self.start) + " " + \
					str(self.end) + " " + str(self.value);
		try:
			return outstring;
		except:
			sys.stderr.write("No BED all information for %s\n" % self)
			return ''


#----------------------------
#----------------------------



class FIXED_STEP:
	"""
	Class to deal with files in fixed step format
	"""
	def __init__(self, position, value=0):
		self.position = position;
		self.value = value;
	def __set__(self, position, value=0):
		self.position = position;
		self.value = value;
	def getAll(self):
		outstring = str(self.position) + " " + str(self.value)
		try:
			return outstring;
		except:
			sys.stderr.write("No information for %s\n" % self)
			return ''


class mem_FIXED_STEP:
	"""
	Class to dead with files in fixed step format, in a memory
	efficient manner
	"""
	def __init__(self, start, step, values):
		self.start = start;
		self.step = step;
		self.values = values;
	def __set__(self, start, step, values):
		self.start = start;
		self.step = step;
		self.values = values;
	def getInfo(self):
		outstring = str(self.start) + " " + str(self.step);




class PAIR:
	"""
	Class to deal with paired-end read bed files
	"""
	def __init__(self, one, two):
		self.one = one;
		self.two = two;
	def __set__(self, one, two):
		self.one = one;
		self.two = two;


#----------------------------
#----------------------------    
    
    
class BED_annotated:
	"""
    Class to deal with bed-like file, assign the needed fields to corresponding variables, and dump the rest into a field called annotation. The bed-like format has to be like BED for the first 6 columns, but there can be more columns after that. 
    """
	def __init__(self, species, file=None, bed_type="BED3", val_threshold= -1):
		""" 
		Overload __init__ so that if a threshold is given, only
		grab bed vals that are above threshold -- This won't do
		anything different for bed3 entries because there is no value
		information in these.
			
		Reads in a bed file and builds a dictionary with chromosomes
		as keys and lists of bed elements as values, 
		Each bed element is a tuple!
		"""
		
		self.bed_vals = {}
		"""
        initialize a dictionary with chromosomes
        """
		for c in GenomeData.species_chroms[species]:
			self.bed_vals[c] = [];

		if(file):
			infile = open(file);
			for line in infile:
				""" check to make sure not a header line """
				if not re.match(comment, line):
					line = line.strip()
					sline = line.split()
					if re.match(bed_type, "BED2"):
						#If want BED2 type: start, strand, chrom is not needed as it is in the dictionary key
						assert (len(sline) >= 6)
						start = int(sline[1])
						strand = sline[5]
						annotation = ""
						if len(sline) >= 7:
							annotation = "\t".join(sline[2:5] + sline[6:])
						if (sline[0] in self.bed_vals.keys()) and (float(sline[4]) >= val_threshold):
							bed = (start, strand, annotation)
							self.bed_vals[sline[0]].append(bed)
					elif re.match(bed_type, "BED3"):#Class for bed lines with 3 values: chrom, start and end
						assert (len(sline) >= 3)
						chrom = sline[0]
						start = int(sline[1])
						end =  int(sline[2])
						annotation = ""
						if len(sline) >= 4:
							annotation = "\t".join(sline[3:])
						if sline[0] in self.bed_vals.keys():
							bed = (chrom, start, end, annotation)
							self.bed_vals[sline[0]].append(bed)
					elif re.match(bed_type, "BED4"):#start, end, name, strand
						assert (len(sline) >= 6)
						#sline[0] is chrom
						start = int(sline[1])
						end =  int(sline[2])
						name = sline[3]
						value = float(sline[4])
						strand = sline[5]
						annotation = ""
						if len(sline) >= 7:
							annotation = "\t".join(sline[4] + sline[6:])
						if (sline[0] in self.bed_vals.keys()) and (value >= val_threshold):
							bed = (start, end, name, strand, annotation)
							self.bed_vals[sline[0]].append(bed)
					elif re.match(bed_type, "BED_GRAPH"):# chrom, start, end, value
						assert (len(sline) >= 5)
						chrom = sline[0]
						start = int(sline[1])
						end =  int(sline[2])
						value =  float(sline[4])
						annotation = ""
						if len(sline) >= 6:
							annotation = "\t".join(sline[3] + sline[5:])
						if (sline[0] in self.bed_vals.keys()) and (value >= val_threshold):
							bed = (chrom, start, end, value, annotation)
							self.bed_vals[sline[0]].append(bed);
					elif re.match(bed_type, "BED6"):#chrom, start, end, name, score, strand
						assert (len(sline) >= 6)
						chrom = sline[0]
						start = int(sline[1])
						end =  int(sline[2])
						name = sline[3]
						value =  float(sline[4])
						strand = sline[5]
						annotation = ""
						if len(sline) >= 7:
							annotation = "\t".join(sline[6:])
						if (sline[0] in self.bed_vals.keys()) and (value >= val_threshold):
							bed = (chrom, start, end, name, value, strand, annotation)
							self.bed_vals[sline[0]].append(bed);
							
#----------------------------
#----------------------------

class BED:
	"""
	Class to deal with bed  and bed_graph files and do common operations 
	Key entity: 
	bed_vals: {chrom:[BED2/BED_GRAPH/... objects]}
	""" 
    
	def __init__(self, species="hg18", file=None, bed_type="BED3", val_threshold = -1):

		""" 
		Overload __init__ so that if a threshold is given, only
		grab bed vals that are above threshold -- This won't do
		anything different for bed3 entries because there is no value
		information in these.
			
		Reads in a bed file and builds a dictionary with chromosomes
		as keys and lists of bed elements as values
		"""
		
		self.bed_vals = {}

		"""
		initialize a dictionary with chromosomes
		"""
		for c in GenomeData.species_chroms[species]:
			self.bed_vals[c] = [];

		if(file):
			if re.match(bed_type, "BED3"):	
				infile = open(file);
				for line in infile:
					""" check to make sure not a header line """
					if not re.match("track", line):
						line = line.strip();
						sline = line.split();
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file
						bed = BED3(sline[0], atoi(sline[1]), atoi(sline[2]));
						if sline[0] in self.bed_vals.keys():
							self.bed_vals[sline[0]].append(bed);
						#if len(sline) == 3:
							#bed = BED3(sline[0], atoi(sline[1]), atoi(sline[2]));
							#self.bed_vals[sline[0]].append(bed);
						#if len(sline) == 4:
							#if atof(sline[3]) >= val_threshold:
								#bed = BED3(sline[0], atoi(sline[1]),
										#atoi(sline[2]));
								#self.bed_vals[sline[0]].append(bed);
						#if len(sline) == 6:
							#if atof(sline[4]) >= val_threshold:
								#bed = BED3(sline[0], atoi(sline[1]),
										#atoi(sline[2]));
								#self.bed_vals[sline[0]].append(bed);
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file

			elif re.match(bed_type, "BED4"):
				infile = open(file);
				for line in infile:
					if not re.match("track", line):
						line = line.strip();
						sline = line.split();
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file
						# start, end, name, strand
						bed = BED4(atoi(sline[1]), atoi(sline[2]), sline[3], sline[5]);
						if sline[0] in self.bed_vals.keys():
							self.bed_vals[sline[0]].append(bed);
						#if len(sline) == 6:
							#bed = BED4(atoi(sline[1]), atoi(sline[2]),
									#sline[3], sline[5])
							#self.bed_vals[sline[0]].append(bed);
											
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file

			elif re.match(bed_type, "BED_GRAPH"):
				""" If want BED_GRAPH type """
				infile = open(file);
				for line in infile:
					""" check to make sure not a header line """
					if not re.match("track", line):
						line = line.strip();
						sline = line.split();
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file
						if sline[0] in self.bed_vals.keys():
							if len(sline) <= 3:
								sys.stderr.write("Can't make bed_graph with less than 4 elements");
								raise bedError;
							elif len(sline) == 4: # This format does not conform to the BED format
								if atof(sline[3]) >= val_threshold:
									bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]), atof(sline[3]));
									self.bed_vals[sline[0]].append(bed);
							else:
								if atof(sline[4]) >= val_threshold:
									bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]), atof(sline[4]));
									self.bed_vals[sline[0]].append(bed);
									#if len(sline) == 3:
										#sys.stderr.write("Can't make bed_graph with only \
										#3 elements")
										#raise bedError 
									#if len(sline) == 4:
										#if atof(sline[3]) >= val_threshold:
											#bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]),
													#atof(sline[3]));
											#self.bed_vals[sline[0]].append(bed);
									#if len(sline) == 6:
										#if atof(sline[4]) >= val_threshold:
											#bed = BED_GRAPH(sline[0], atoi(sline[1]), atoi(sline[2]),
													#atof(sline[4]));
											#self.bed_vals[sline[0]].append(bed);
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file


			elif re.match(bed_type, "BED2"):
				""" If want BED2 type """
				infile = open(file);
				for line in infile:
					""" check to make sure not a header line """
					if not re.match("track", line):
						line = line.strip();
						sline = line.split();
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file
						if sline[0] in self.bed_vals.keys():
							if len(sline) <= 5:
								sys.stderr.write("Need BED6 to make BED2")
								raise bedError
							else: # start depends on orientation
								if atof(sline[4]) >= val_threshold:
									if plus.match(sline[5]):
										bed = BED2(atoi(sline[1]), sline[5]);
									elif minus.match(sline[5]):# The BED format is open ended, hence the real start is atoi(sline[2]) - 1
										bed = BED2(atoi(sline[2])-1, sline[5]);
									self.bed_vals[sline[0]].append(bed);
									#if len(sline) == 3:
										#sys.stderr.write("Need BED6 to make BED2")
										#raise bedError
									#if len(sline) == 4:
										#sys.stderr.write("Need BED6 to make BED2")
										#raise bedError
									#if len(sline) == 6:
										#if atof(sline[4]) >= val_threshold:
											#if plus.match(sline[5]):
												#bed = BED2(atoi(sline[1]), sline[5]);
											#elif minus.match(sline[5]):
												#bed = BED2(atoi(sline[2]), sline[5]);
											#self.bed_vals[sline[0]].append(bed);
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file

			elif re.match(bed_type, "BED6"):
				""" If want BED6 element """
				infile = open(file);
				for line in infile:
					""" check to make sure not a header line """
					if not re.match("track", line):
						line = line.strip();
						sline = line.split();
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file
						if sline[0] in self.bed_vals.keys():
							if len(sline) <= 5:
								sys.stderr.write("Need at least 6 columns to make BED6")
								raise bedError
							else:
								if atof(sline[4]) >= val_threshold:
									bed = BED6(sline[0], atoi(sline[1]), atoi(sline[2]), sline[3], atof(sline[4]), sline[5]);
									self.bed_vals[sline[0]].append(bed);
									#if len(sline) == 3:
										#sys.stderr.write("Need BED6 to make BED6")
										#raise bedError
									#if len(sline) == 4:
										#sys.stderr.write("Need BED6 to make BED6")
										#raise bedError
									#if len(sline) == 6:
										#if atof(sline[4]) >= val_threshold:
											#bed = BED6(sline[0], atoi(sline[1]), atoi(sline[2]),
													#sline[3], atof(sline[4]), sline[5]);
											#self.bed_vals[sline[0]].append(bed);
#==== Changed by Weiqun 4/2010. I am assuming the right ordering of columns in the bed file

			# formating might be off, might not work, 11/1/2011
			elif re.match(bed_type, "FIXED_STEP"):
				#If want FIXED_STEP type
				infile = open(file);
				for line in infile:
					line = line.strip();
					sline = line.split();
					##if re.match(sline[0], "fixedStep"):
					## TEST: I think this will be faster than the re.match
					## on every line 
					if sline[0][0] == "f":
						chrom = sline[1].split("=")[1];
						start = atoi(sline[2].split("=")[1]);
						step_size = atoi(sline[3].split("=")[1]);
						position = start;
					else:
						bed = FIXED_STEP(position, atof(sline[0]));
						position += step_size;
						self.bed_vals[chrom].append(bed);

			# formating might be off, might not work, 11/1/2011
			elif re.match(bed_type, "mem_FIXED_STEP"):
				""" fixed step file """
				infile = open(file);
		
				line = infile.readline();
				line = line.strip();
				sline = line.split();
						
				chrom = sline[1].split("=")[1];
				start = atoi(sline[2].split("=")[1]);
				step = atoi(sline[3].split("=")[1]);
				values = [];
		
				for line in infile:
					line = line.strip();
					sline = line.split();
		
					##if re.match(sline[0], "fixedStep"):
					## TEST: I think this will be faster than the re.match
					## on every line 
					if sline[0][0] == "f":
						bed = mem_FIXED_STEP(start, step, values);
						self.bed_vals[chrom].append(bed);
						chrom = sline[1].split("=")[1];
						start = atoi(sline[2].split("=")[1]);
						step = atoi(sline[3].split("=")[1]);
						values = [];
					else:
						values.append(atof(sline[0]));
						bed = mem_FIXED_STEP(start, step, values);
						self.bed_vals[chrom].append(bed);


    
	def keys(self):
		"""
		Return a list of the keys - duplicating the function of a dictionary
		"""
		return self.bed_vals.keys()

	"""
	-- BED2, BED3, BED_GRAPH, FIXED_STEP all ok -- BED3, BED_GRAPH AND
	FIXED STEP DO NOT HAVE STRAND INFO SO DON'T NEEED TO WORRY ABOUT
	THEM HERE
	-- BED2 STRANDS WERE TAKE INTO CONSIDERATION WHEN READ IN
	-- BED6 - NEED TO USE getStarts_consider_strands
	"""

	def getStarts_consider_strands(self, chrom):
		"""
		Return a list of starts on a given chromosome
		"""
		starts = [];
		for t in self.bed_vals[chrom]:
			if plus.match(t.strand):
				starts.append(t.start);
			elif minus.match(t.strand):
				starts.append(t.end);
		try:
			return starts;
		except:
			sys.stderr.write("Having trouble returning starts %s\n" % self)
			return ''

	def getStarts(self, chrom):
		"""
		Return a list of starts on a given chromosome
		"""
		starts = [];
		for t in self.bed_vals[chrom]:
			starts.append(t.start);
		try:
			return starts;
		except:
			sys.stderr.write("Having trouble returning starts %s\n" % self)
			return ''


	def getEnds(self, chrom):
		"""
		Return a list of starts on a given chromosome
		"""
		ends = [];
		for t in self.bed_vals[chrom]:
			ends.append(t.end);
		try:
			return ends;
		except:
			sys.stderr.write("Having trouble returning ends %s\n" % self)
			return ''

	def getChroms(self):
		"""
		Return a list of all the chromosomes in order (ordered keys)
		"""
		try:
			return chroms[species];
		except:
			sys.stderr.write("Can't return chromosome list\n" % self)
			return ''


	def getNumVals(self):
		"""
		Return the number of bed vals in BED instance
		"""
		num = 0;
		for c in self.bed_vals.keys():
			num += len(self.bed_vals[c]);
		return num;


	def addChrom(self, chrom, bed_list):
		if self.bed_vals.has_key(chrom) and len(self.bed_vals[chrom]) > 0:
			sys.stderr.write("chromsome %s already populated\n" % chrom)
			raise bedError
		else: self.bed_vals[chrom] = bed_list;
		

	def __del__(self):
		"""
		Delete, delete;
		"""
		self.bed_vals.clear()


	def __contains__(self, item):
		"""
		Returns  mapping iterator
		"""
		return self.bed_vals.has_key(item)

	def __iter__(self):
		"""
		Returns mapping iterator
		"""
		return self.bed_vals.iterkeys()

	def __len__(self):
		"""
		Returns number of bed_vals
		"""
		return len(self.bed_vals)

	def __delitem__(self, name):
		"""
		removes a chrom if its name exists in the dictionary
		-- I guess this could possible be useful at some point
		
		"""
		if self.bed_vals.has_key(name):
			del self.bed_vals[name]
		else: raise bedError
		
	def __setitem__(self, name, bedlist):
		"""
		Sets a new bed value
		"""
		self.bed_vals[name] = bedlist

	def __getitem__(self, name):
		"""
		Returns a bed_val indexed by its name or None if no such bed_val exists
		"""
		if self.bed_vals.has_key(name):
			return self.bed_vals[name]
		else: raise bedError



class BED_PAIR:
	"""
	Class to deal with bed files and do common operations 
	""" 
	def __init__(self, species="hg18", file=None, bed_type="BED2", val_threshold=-1e308):

		""" 
		Overload __init__ so that if a threshold is given, only
		grab bed vals that are above threshold -- This won't do
		anything different for entries with no value information
			
		Reads in a bed file and builds a dictionary with chromosomes
		as keys and lists of bed elements as values

		** RIGHT NOW ONLY DOING BED2 **
		"""
		
		self.bed_pairs = {}
		"""
		initialize a dictionary with chromosomes
		"""
		for c in GenomeData.species_chroms[species]:
			self.bed_pairs[c] = [];
		if(file):
			if re.match(bed_type, "BED2"):
				""" 
				store pair information - right now only storing
				BED2 information for each individual read
				
				IMPORTANT: This assumes that pairs have aleady been filtered
				to only keep those pairs where:
				(1) both are on the same chromosome
				(2) they are properly oriented (i.e. facing each other)
				"""

				infile = open(file);
				for line in infile:
					if not re.match("track", line):
						line = line.strip();
						pline = line.split("PR");

						chrom = line.split()[0];
						r_one = pline[0].split();
						r_two = pline[1].split();

						if plus.match(r_one[5]):
							bed_one = BED2(atoi(r_one[1]), r_one[5]);
							bed_two = BED2(atoi(r_two[2]), r_two[5]);
							bed_pair = PAIR(bed_one, bed_two);
						
						elif minus.match(r_one[5]):
							bed_one = BED2(atoi(r_one[2]), r_one[5])
							bed_two = BED2(atoi(r_two[1]), r_two[5]);
							bed_pair = PAIR(bed_two, bed_one);

						self.bed_pairs[chrom].append(bed_pair);


	def keys(self):
		"""
		Return a list of the keys - duplicating the function of a dictionary
		"""
		return self.bed_pairs.keys()


	def getChroms(self):
		"""
		Return a list of all the chromosomes in order (ordered keys)
		"""
		try:
			return chroms[species];
		except:
			sys.stderr.write("Can't return chromosome list\n" % self)
			return ''


	def getNumVals(self):
		"""
		Return the number of bed vals in BED instance
		"""
		num = 0;
		for c in self.bed_pairs.keys():
			num += len(self.bed_pairs[c]);
		return num;


	def addChrom(self, chrom, bed_list):
		if self.bed_pairs.has_key(chrom) and len(self.bed_pairs[chrom]) > 0:
			sys.stderr.write("chromsome %s already populated\n" % chrom)
			raise bedError
		else: self.bed_pairs[chrom] = bed_list;
		

	def __del__(self):
		"""
		Delete, delete;
		"""
		self.bed_pairs.clear()


	def __contains__(self, item):
		"""
		Returns  mapping iterator
		"""
		return self.bed_pairs.has_key(item)

	def __iter__(self):
		"""
		Returns mapping iterator
		"""
		return self.bed_pairs.iterkeys()

	def __len__(self):
		"""
		Returns number of bed_pairs
		"""
		return len(self.bed_pairs)

	def __delitem__(self, name):
		"""
		removes a chrom if its name exists in the dictionary -- I
		guess this could possible be useful at some point
		"""
		if self.bed_pairs.has_key(name):
			del self.bed_pairs[name]
		else: raise bedError
		
	def __setitem__(self, name, bedlist):
		"""
		Sets a new bed value
		"""
		self.bed_pairs[name] = bedlist

	def __getitem__(self, name):
		""" Returns a bed_pair indexed by its name or None if no such
		bed_pair exists """
		if self.bed_pairs.has_key(name):
			return self.bed_pairs[name]
		else: raise bedError
