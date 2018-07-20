#!/usr/bin/env python
# Authors: Dustin E Schones, Weiqun Peng and Keji Zhao
#
# Extended by Weiqun Peng 2009.7, to add TSS, genebody, promoter+gene body functionality
#
# 2011.5
# modified getGenebodys(self, downstream) so that genebodies whose length are shorted than downstream are ignored. 
# added : 	def getClosestTSS(self, chrom, position):	
		#def getClosestTES(self, chrom, position):
		#def getIntergenicRegion(self, upstreamfar, upstreamclose):

# 2011.9
# Added exon function


# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# schonesde@mail.nih.gov)
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

"""
This module contains the classes to deal with UCSC
known gene files
"""
import re, os, shutil, time, sys
from math import *
from string import *
import bisect
from operator import itemgetter

plus = re.compile('\+');
minus = re.compile('\-');

import BED;
import Utility_extended
UCSCError = "Error in UCSC class";


class UCSC:
	"""
	Class for keeping known transcript information.  
	additional annotations: gene symbol, Entriz ID, etc
	internalID: a unique identifier for each entry. It is defined inside. 
	"""

	def __init__(self, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, additional_annotations=[], internalID=0):
		self.name = name;
		self.chrom = chrom;
		self.strand = strand;
		self.txStart = txStart;
		self.txEnd = txEnd;
		self.cdsStart = cdsStart;
		self.cdsEnd = cdsEnd;
		self.exonCount = exonCount;
		self.exonStarts = exonStarts;
		self.exonEnds = exonEnds;
		self.additional_annotations = additional_annotations
		self.internalID = internalID
		
	def __set__(self, name, chrom, strand, txStart, txEnd, cdsStart,
				cdsEnd, exonCount, exonStarts, exonEnds, additional_annotations=[], internalID=0):
		self.name = name;
		self.chrom = chrom;
		self.strand = strand;
		self.txStart = txStart;
		self.txEnd = txEnd;
		self.cdsStart = cdsStart;
		self.cdsEnd = cdsEnd;
		self.exonCount = exonCount;
		self.exonStarts = exonStarts;
		self.exonEnds = exonEnds;
		self.additional_annotations = additional_annotations
		self.internalID = internalID

	def getAll(self):
		outstring = self.name + "\t" + self.chrom + "\t" + self.strand + "\t" + \
				str(self.txStart) + \
					"\t" + str(self.txEnd) + "\t" + str(self.cdsStart) + "\t" + \
					str(self.cdsEnd) + "\t" + str(self.exonCount) + "\t" + \
					str(self.exonStarts) + "\t" + str(self.exonEnds) + "\t" + \
					"\t".join([str(i) for i in self.additional_annotations]) + "\n"
		try:
			return outstring;
		except:
			sys.stderr.write("No UCSC known gene information for %s\n" % self)
			return ''
	
	def getPromoter(self, upstream, downstream):
		"""
		Return the promoter coordinates defined by given upstream and
		downstream distances from TSS (txStart) 
		
		If a transcript is so short that the promoter end is beyond TES, this transcript will be discarded.
	
		NOTICE:  if transcript is on negative strand, the promoter start
		and end are still sequential
		"""
		include = 1
		if plus.match(self.strand):
			prom_start = self.txStart - upstream;
			prom_end = self.txStart + downstream;
			if prom_start < 0:
				include = 0
			if prom_end >= self.txEnd:
				include = 0
		elif minus.match(self.strand):
			prom_start = self.txEnd - downstream;
			prom_end = self.txEnd + upstream;
			if prom_start <= self.txStart:
				include = 0
		if include > 0:
			return (prom_start, prom_end);
		else:
			return ();
		
	def get5UTR(self, upstream, downstream):
		"""
		Return the 5UTR coordinates defined by given upstream and
		downstream distances from the boudaries 
		
		If a transcript is so short that the 5UTR end is beyond TES, this transcript will be discarded.
	
		NOTICE:  if transcript is on negative strand, the 5UTR start
		and end are still sequential
		"""
		include = 1
		if plus.match(self.strand):
			fUTR_start = self.txStart - upstream;
			fUTR_end = self.cdsStart + downstream;
			if fUTR_start < 0 or fUTR_end >= self.txEnd:
				include = 0
		elif minus.match(self.strand):
			fUTR_start = self.cdsEnd - downstream;
			fUTR_end = self.txEnd + upstream;
			if fUTR_start <= self.txStart:
				include = 0
		if include > 0:
			return (fUTR_start, fUTR_end)
		else:
			return ()

	def getExons(self):
		"""
		These are internal exons, the 3'UTR and 5'UTR do not count
		exons are represented by a list of [(start,end)], start <= end, the output shout be self sorted already
		"""
		exons=[]
		if self.exonCount > 0:
			exon_Starts_str = (self.exonStarts.split(','))[:-1] #remove the last '' because the format is '1,2,3,'
			exon_Ends_str = (self.exonEnds.split(','))[:-1]#remove the last '' because the format is '1,2,3,'
			exon_Starts = [int(x) for x in exon_Starts_str]
			exon_Ends = [int(x) for x in exon_Ends_str]
			assert len(exon_Starts) == len(exon_Ends)
			exons = zip(exon_Starts, exon_Ends)
		return exons;

	def getIntrons(self):
		"""
		introns are represented by a list of [(start,end)], start <= end, the output should be self-sorted already
		"""
		introns=[]
		if self.exonCount > 1:
			exon_Starts_str = (self.exonStarts.split(','))[:-1] #remove the last '' because the format is '1,2,3,'
			exon_Ends_str = (self.exonEnds.split(','))[:-1]#remove the last '' because the format is '1,2,3,'
			exon_Starts = [int(x) for x in exon_Starts_str]
			exon_Ends = [int(x) for x in exon_Ends_str]
			intron_Starts = [x+1 for x in exon_Ends[:len(exon_Ends)-1]]
			intron_Ends = [x-1 for x in exon_Starts[1:]]
			assert len(intron_Ends) == len(intron_Starts)
			for i in xrange(len(intron_Ends)):
				assert intron_Starts[i] <= intron_Ends[i]
			introns = zip(intron_Starts, intron_Ends)
		return introns

	def get3UTR(self, upstream, downstream):
		"""
		Return the 3UTR coordinates defined by given upstream and
		downstream distances from the boudaries (start, end)
		
		If a transcript is so short that the 3UTR end is beyond TSS, this transcript will be discarded.
	
		NOTICE:  if transcript is on negative strand, the 3UTR start
		and end are still sequential
		"""
		include = 1
		if plus.match(self.strand):
			tUTR_start = self.cdsEnd - upstream;
			tUTR_end = self.txEnd + downstream;
			if tUTR_start < self.txStart:
				include = 0
		elif minus.match(self.strand):
			tUTR_start = self.txStart - downstream;
			tUTR_end = self.cdsStart + upstream;
			if tUTR_start < 0 or tUTR_end > self.txEnd:
				include = 0
		if include > 0:
			return (tUTR_start, tUTR_end)
		else:
			return ()
		
	def getGeneEnd(self, upstream, downstream):
		"""
		Return the GeneEnd coordinates defined by given upstream and
		downstream distances from TES (txEnd) 
		
		If a transcript is so short that the promoter end is beyond TSS, this transcript will be discarded.
	
		NOTICE:  if transcript is on negative strand, the promoter start
		and end are still sequential
		"""
		include = 1
		if plus.match(self.strand):
			ge_start = self.txEnd - upstream;
			ge_end = self.txEnd + downstream;
			if ge_start < self.txStart:
				include = 0
		elif minus.match(self.strand):
			ge_start = self.txStart - downstream;
			ge_end = self.txStart + upstream;
			if ge_start < 0 or ge_end > self.txEnd:
				include = 0
		if include > 0:
			return (ge_start, ge_end);
		else:
			return ()
	
	def getExtendedGeneBody(self, upstream, downstream):
		"""
		returns a tuple (egb_start, egb_end)
		"""
		if plus.match(self.strand):
			egb_start = self.txStart - upstream;
			egb_end = self.txEnd + downstream;
		elif minus.match(self.strand):
			egb_start = self.txStart - downstream;
			egb_end = self.txEnd + upstream;
		return (egb_start, egb_end);
	
	def Is3UTRIntronContaining(self, allowance=0):
		"""
		allowance is for what?
		"""
		yes = 0
		introns = self.getIntrons() #[(start, end)]
		if len(introns) > 0:
			if plus.match(self.strand):
				start = self.cdsEnd
				end = self.txEnd
				last_intron = introns[-1]
			elif minus.match(self.strand):
				start = self.txStart
				end = self.cdsStart
				last_intron = introns[0]
			start = min(start + allowance, end)
			end = max(start, end - allowance)
			if Utility_extended.overlap(start, end, last_intron[0], last_intron[1]) == 1: 
				yes = 1
		return yes

#----------------------------
#----------------------------

class UCSC_lite:   
    """
    Class for keeping known gene information.  Only partial (most
    usefull) information stored.
    """
    def __init__(self, name, chrom, strand, txStart, txEnd):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;
    def __set__(self, name, chrom, strand, txStart, txEnd):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;

    def getAll(self):
        outstring = self.name + "\t" + self.chrom + "\t" + self.strand + \
                    "\t" + str(self.txStart) + "\t" + str(self.txEnd);
        try:
            return outstring;
        except:
            sys.stderr.write("No UCSC known gene information for %s\n" % self)
            return ''


#----------------------------
#----------------------------


class KnownGenes:
	"""
	Class to read in UCSC known gene files and get all the info
	"""
	
	def __init__(self, file=None):
		"""
		Reads in a UCSC known gene file and builds a dictionary
		with chromosomes as keys and lists representing
		the transcripts on each chromosome as values.
		
		"""
		self.gene_coords = {}
		infile = open(file);
		index = 0
		for line in infile:
			""" check to make sure not a header line """
			if not re.match("#", line):
				line = line.strip();
				sline = line.split("\t");
				if sline[1] not in self.gene_coords.keys(): #sline[1] is chrom
					self.gene_coords[sline[1]] = [];
				if len(sline) == 10:
					coord = UCSC(sline[0], sline[1], sline[2], atoi(sline[3]), atoi(sline[4]), atoi(sline[5]), atoi(sline[6]), atoi(sline[7]), sline[8], sline[9], [], index);
				else:
					coord = UCSC(sline[0], sline[1], sline[2], atoi(sline[3]), atoi(sline[4]), atoi(sline[5]), atoi(sline[6]), atoi(sline[7]), sline[8], sline[9], sline[10:], index);
				self.gene_coords[sline[1]].append(coord);
				index += 1

	def getPromoters(self, upstream, downstream):
		"""
		If a gene is so short that the promoter end is beyond TES, this gene will be discarded.
	
		NOTICE:  even if gene is on negative strand, the promoter start
		and end are still sequential
		
		Return:
		{chrom:[UCSC_lite]}
		UCSC_lite: name, chrom, strand, start, end
		"""
		self.prom_coords = {};
		include = 1
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				cord = g.getPromoter(upstream, downstream)
				if len(cord) > 0:
					ucsc = UCSC_lite(g.name, chrom, g.strand, cord[0], cord[1]);
					if chrom not in self.prom_coords.keys():
						self.prom_coords[chrom] = [];
					self.prom_coords[chrom].append(ucsc);
		return self.prom_coords;

	def IsInPromoter(self, chrom, position, strand, upstream, downstream):
		"""
		strand can be '+', '-', "ALL"
		
		If the position is not strand specific:
			returns the list of transcripts (ucsc objects) whose promoters encompass the position
		If the position is strand specific: 
			If the position is located in a promoter in the same strand, returns the transcripts (ucsc objects)
		The reason for returning the ucsc objects is that the same gene name can be used for multiple entries.
		"""
		transcripts = []
		if len(self.gene_coords[chrom]) > 0:
			if strand.strip() == 'ALL':
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart - upstream
						end = g.txStart + downstream
						if (position >= max(start,0)) and (position <= min(end, g.txEnd)):
							transcripts.append(g)
					elif minus.match(g.strand):
						start = g.txEnd - downstream
						end = g.txEnd + upstream
						if position >= max(start, g.txStart) and position <= end:
							transcripts.append(g)
					else: 
						print "Gene strand identification crisis!";
						exit(1)	
			elif plus.match(strand):
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart - upstream
						end = g.txStart + downstream
						if (position >= max(start,0)) and (position <= min(end, g.txEnd)):
							transcripts.append(g)
			elif minus.match(strand):
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						start = g.txEnd - downstream
						end = g.txEnd + upstream
						if position >= max(start, g.txStart) and position <= end:
							transcripts.append(g)
		return transcripts

	def getTSS(self):
		"""
		Return a dictionary of TSS positions keyed by gene name.  
		"""
		self.TSS={};
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if plus.match(g.strand):
					TSS[g.name]=g.txStart
				elif minus.match(g.strand):
					TSS[g.name]=g.txEnd
		return self.TSS
				
				
				
	def getClosestTSS(self, chrom, position, strand):
		"""
		Find the gene whose TSS is the closest to the specified position. 
		If there is degeneracy, pick the left one?
		position needs to be within the boundary of chromosome, checking of that is done outside.
		returns the gene_name and distance
		"""
		if len(self.gene_coords[chrom]) > 0:
			if plus.match(strand):
				g = (self.gene_coords[chrom])[0]
				if plus.match(g.strand):
					mymin = abs(g.txStart - position)
					gene_name = g.name 
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						if mymin > abs(g.txStart - position):
							mymin = abs(g.txStart - position)
							gene_name = g.name
			elif minus.match(strand):
				if minus.match(g.strand):
					mymin = abs(g.txEnd - position)
					gene_name = g.name	
				for g in self.gene_coords[chrom]:	
					if minus.match(g.strand):
						if mymin > abs(g.txEnd - position):
							mymin = abs(g.txEnd - position)
							gene_name = g.name	
			else: print "Gene strand identification crisis!";	
			return (gene_name, mymin)
		else: 
			print "There is no gene on ", chrom;
			return ()

	def get5UTRs(self, upstream, downstream):
		"""
		Return the 5UTR coordinates defined by given upstream and
		downstream distances from the boundaries in a bed dictionary
		
		If a gene is so short that the 5UTR end is beyond TES, this gene will be discarded.
	
		NOTICE:  even if gene is on negative strand, the 5UTR start
		and end are still sequential
		"""
		self.fUTR_coords = {};
		include = 1
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				cord = g.get5UTR(upstream, downstream)
				if len(cord) > 0:
					ucsc = UCSC_lite(g.name, chrom, g.strand, cord[0], cord[1]);
					if chrom not in self.fUTR_coords.keys():
						self.fUTR_coords[chrom] = [];
					self.fUTR_coords[chrom].append(ucsc);
		return self.fUTR_coords;


	def getGenebodys(self, downstream):
		"""
		downstream: downstream of TSS
		Gene body and promoter are mutually exclusive
		Returns: {chrom:[UCSC_lite]} 
		
		NOTICE:  if gene is on negative strand, the region start
		and end are still sequential
		"""
		self.genebody={};
		for chrom in self.gene_coords.keys():
			self.genebody[chrom]=[]
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if plus.match(g.strand):
					GB_start = g.txStart + downstream
					GB_end = g.txEnd
					if GB_start < g.txEnd: # if the gene body is shorter than downstream, ignore
						#print g.name, g.txStart, g.txEnd;
						ucsc = UCSC_lite(g.name, chrom, g.strand, GB_start, GB_end);
						self.genebody[chrom].append(ucsc);
				# When strand is negative, the regions are still sequential
				elif minus.match(g.strand):
					GB_start = g.txStart
					GB_end = g.txEnd - downstream
					if GB_end > GB_start: # if the gene body is shorter than downstream, ignore
						#print g.name, g.txStart, g.txEnd;
						ucsc = UCSC_lite(g.name, chrom, g.strand, GB_start, GB_end);
						self.genebody[chrom].append(ucsc);
				else:
					print "g.strand is", g.strand
					print "Gene strand identification crisis!";
		return self.genebody;

	def IsInGenebody(self, chrom, position, strand, downstream):
		"""
		strand can be '+', '-', "ALL"
		If the position is located in a genebody, returns the gene 
		"""
		transcripts = []
		
		if strand.strip() == 'ALL':
			if len(self.gene_coords[chrom]) > 0:
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart + downstream
						end = g.txEnd
						if position >= start and position <= end:
							transcripts.append(g)
					elif minus.match(g.strand):
						start = g.txStart
						end = g.txEnd - downstream
						if position >= start and position <= end:
							transcripts.append(g)
					else: 
						print "Gene strand identification crisis!";
						exit(1)	
			elif plus.match(strand):
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart + downstream
						end = g.txEnd
						if (position >= max(start,0)) and (position <= end):
							transcripts.append(g)
			elif minus.match(strand):
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						start = g.txStart
						end = g.txEnd - downstream
						if position >= start and position <= end:
							transcripts.append(g)
			else:
				print "Warning: Variable strand wrong input." 
				
		return 	transcripts			

	def getPromotergenebodys(self, upstream):
		"""
		Gene body and promoter are mutually exclusive
		Returns:
			{chrom:[UCSC_lites]}
		
		NOTICE:  if gene is on negative strand, the region start and end are still sequential
		"""
		self.pg={};
		for chrom in self.gene_coords.keys():
			self.pg[chrom]=[];
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if plus.match(g.strand):
					pg_start = g.txStart - upstream
					if pg_start <0:
						pg_start = 0;
					pg_end = g.txEnd
				# When strand is negative, the regions are still sequential
				elif minus.match(g.strand):
					pg_start = g.txStart
					pg_end = g.txEnd + upstream
				else: print "Gene strand identification crisis!";
				ucsc = UCSC_lite(g.name, chrom, g.strand, pg_start, pg_end);
				self.pg[chrom].append(ucsc);
		return self.pg;


	def getExtendedGeneBodys(self, upstream, downstream):
		"""
		Returns a dictionary keyed by chrom, each chrom has a list of ucsc_lite objects.
		NOTICE:  even if gene is on negative strand, the promoter start
		and end are still sequential
		"""
		self.egb_coords = {};
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				cord = g.getExtendedGeneBody(upstream, downstream)
				if len(cord) > 0:
					ucsc = UCSC_lite(g.name, chrom, g.strand, cord[0], cord[1]);
					if chrom not in self.egb_coords.keys():
						self.egb_coords[chrom] = [];
					self.egb_coords[chrom].append(ucsc);
		return self.egb_coords

	def IsInPromoterGenebody(self, chrom, position, strand, upstream):
		"""
		If the position is located in a promoter+genebody, returns the gene name
		"""
		transcripts = []
		if len(self.gene_coords[chrom]) > 0:
			if strand.strip() == 'ALL':
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart - upstream
						end = g.txEnd
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
					elif minus.match(g.strand):
						start = g.txStart
						end = g.txEnd + upstream
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
					else: 
						print "Gene strand identification crisis!";
						exit(1)	
			elif plus.match(strand):
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.txStart - upstream
						end = g.txEnd
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
			elif minus.match(strand):
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						start = g.txStart
						end = g.txEnd + upstream
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
		return transcripts	

	def getExons(self):
		"""
		Return a dictionary (keyed by chrom, values are lists of UCSC_lite objects)
		no matter whether the strand is + or -, the exon start is alway smaller than the exon end 
		"""
		self.exons={}; #{chrom:ucsc_lite}
		for chrom in self.gene_coords.keys():
			self.exons[chrom]=[];
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				the_exons = g.getExons()
				for i in xrange(len(the_exons)):
					exon_name = g.name + '_exon' + str(i);
					ucsc = UCSC_lite(exon_name, chrom, g.strand, the_exons[i][0], the_exons[i][1])
					self.exons[chrom].append(ucsc)
		return self.exons;
	
	def getIntrons(self):
		"""
		Return a dictionary (keyed by chrom, values are lists of UCSC_lite objects)
		no matter whether the strand is + or -, the exon start is alway smaller than the exon end 
		"""
		self.introns={};
		for chrom in self.gene_coords.keys():
			self.introns[chrom]=[];
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				the_introns = g.getIntrons()
				if len(the_introns) > 0:
					for i in xrange(len(the_introns)):
						intron_name = g.name + '_intron' + str(i);
						ucsc = UCSC_lite(intron_name, chrom, g.strand, the_introns[i][0], the_introns[i][1])
						self.introns[chrom].append(ucsc)
		return self.introns;

	def get3UTRs(self, upstream, downstream):
		"""
		Return a dictionary (keyed by chrom, values are lists of UCSC_lite objects)
		"""
		self.threeUTRs={};
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				my3UTR=g.get3UTR(upstream, downstream)
				if len(my3UTR) > 0:
					ucsc = UCSC_lite(g.name, chrom, g.strand, my3UTR[0], my3UTR[1]);
					if chrom not in self.threeUTRs.keys():
						self.threeUTRs[chrom] = [];
					self.threeUTRs[chrom].append(ucsc);
		return self.threeUTRs;


	def IsIn3UTR(self, chrom, position, strand, downstream):
		"""
		The orientation can be '+', '-', and "ALL", the choice of "ALL" is for the peaks without orientations.
		a peak is characterized by its position and its orientation. Example of peaks with orientation is PA-Seq peaks. 
		
		On positive strand: 3UTR = [self.cdsEnd, self.txEnd] 
		On negative strand: 3UTR = [self.txStart, self.cdsStart]
		
		"""
		transcripts = []
		if len(self.gene_coords[chrom]) > 0:
			if strand.strip() == 'ALL':
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.cdsEnd
						end = g.txEnd + downstream
						if position >= start and position <= end:
							transcripts.append(g)
					elif minus.match(g.strand):
						start = g.txStart - downstream
						end = g.cdsStart
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
					else: 
						print "Gene strand identification crisis!";
						exit(1)	
			elif plus.match(strand):
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						start = g.cdsEnd
						end = g.txEnd + downstream
						if position >= start and position <= end:
							transcripts.append(g)
			elif minus.match(strand):
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						start = g.txStart - downstream
						end = g.cdsStart
						if position >= max(start,0) and position <= end:
							transcripts.append(g)
		return transcripts	

	def getTES(self):
		"""
		Return a dictionary of TES positions keyed by gene name.  
		"""
		self.TES={};
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if plus.match(g.strand):
					self.TES[g.name]=g.txEnd
				elif minus.match(g.strand):
					self.TES[g.name]=g.txStart
		return self.TES
		

	def getClosestTES(self, chrom, position, strand):
		"""
		Find the gene whose TES is the closest to the specified position. 
		position needs to be within the boudary of chromosome, checking of that is done outside.
		
		"""
		if len(self.gene_coords[chrom]) > 0:
			if plus.match(strand):
				g = (self.gene_coords[chrom])[0]
				if plus.match(g.strand):
					mymin = abs(g.txEnd - position)
					gene_name = g.name 
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						if mymin > abs(g.txEnd - position):
							mymin = abs(g.txEnd - position)
							gene_name = g.name
			elif minus.match(strand):
				if minus.match(g.strand):
					mymin = abs(g.txStart - position)
					gene_name = g.name	
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						if mymin > abs(g.txStart - position):
							mymin = abs(g.txStart - position)
							gene_name = g.name	
					else: print "Gene strand identification crisis!";	
			else: 
				print "Gene strand identification crisis!";		
			return (gene_name, mymin)
		else: 
			print "There is no gene on ", chrom;
			return ()


	def getGeneEnds(self, upstream, downstream):
		"""
		Return the GeneEnd coordinates defined by given upstream and
		downstream distances from TES (txEnd) in a bed dictionary
		{chrom:[UCSC_lite]}
		
		If a gene is so short that the GeneEnd is beyond TSS, this gene will be discarded.
	
		NOTICE:  even if gene is on negative strand, the GeneEnd start
		and end are still sequential
		"""
		ge_coords = {}
		for chrom in self.gene_coords.keys():
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				cord = g.getGeneEnd(upstream, downstream)
				if len(cord)>1:
					ucsc = UCSC_lite(g.name, chrom, g.strand, cord[0], cord[1]);
					if chrom not in ge_coords.keys():
						ge_coords[chrom] = [];
					ge_coords[chrom].append(ucsc);
		return ge_coords;



	def overlap(self, start1, end1, start2, end2):
		assert start1 <= end1
		assert start2 <= end2
		if end1 < start2 or end2 < start1: #Non-overlap
			return 0
		else:
			return 1 
	
	def IsOverlappingPromoterGenebody(self, chrom, start, end, strand, upstream):
		"""
		If the region overlapps with a promoter+genebody, returns gene
		"""
		transcripts = []
		if len(self.gene_coords[chrom]) > 0:
			if strand.strip() == 'ALL':
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						pg_start = g.txStart - upstream
						pg_end = g.txEnd
						if self.overlap(pg_start, pg_end, start, end) > 0:
							transcripts.append(g)
					elif minus.match(g.strand):
						pg_start = g.txStart
						pg_end = g.txEnd + upstream
						if self.overlap(pg_start, pg_end, start, end) > 0:
							transcripts.append(g)
					else: 
						print "Gene strand identification crisis!";
						exit(1)
			elif plus.match(strand):
				for g in self.gene_coords[chrom]:
					if plus.match(g.strand):
						pg_start = g.txStart - upstream
						pg_end = g.txEnd
						if self.overlap(pg_start, pg_end, start, end) > 0:
							transcripts.append(g)
			elif minus.match(strand):
				for g in self.gene_coords[chrom]:
					if minus.match(g.strand):
						pg_start = g.txStart
						pg_end = g.txEnd + upstream
						if self.overlap(pg_start, pg_end, start, end) > 0:
							transcripts.append(g)
		return transcripts
	
	def getIntergenicRegion(self, upstreamfar, upstreamclose):
		"""
		If the prescribed region overlaps with neighbouring gene, discard it
		Returns a dictionary
		"""
		assert upstreamfar > upstreamclose
		ir={}
		for chrom in self.gene_coords.keys():
			ir[chrom]=[];
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if plus.match(g.strand):
					start = g.txStart - upstreamfar
					end = g.txStart - upstreamclose
					# If no overlap with other transcripts, register
					if start > 0 and len(self.IsOverlappingPromoterGenebody(chrom, start, end, upstreamclose)) == 0 :
						ucsc = UCSC_lite(g.name, chrom, g.strand, start, end);
						ir[chrom].append(ucsc)
				# When strand is negative, the regions are still sequential
				elif minus.match(g.strand):
					start = g.txEnd + upstreamclose
					end = g.txEnd + upstreamfar
					if len(self.IsOverlappingPromoterGenebody(chrom, start, end, upstreamclose)) == 0 :
						ucsc = UCSC_lite(g.name, chrom, g.strand, start, end);
						ir[chrom].append(ucsc)
				else: 
					print "Gene strand identification crisis!";
					exit(1)
		return ir
				
	def lookupByName(self, chrom, name):
		"""
		might return a list of gene
		"""
		transcripts = []
		chrom_coords = self.gene_coords[chrom]
		for gene in chrom_coords:
			if name == gene.name:
				transcripts.append(gene)
		return transcripts
	
	def lookupByInternalID(self, chrom, ID):
		"""
		Should return a unique gene
		"""
		transcripts = []
		chrom_coords = self.gene_coords[chrom]
		for gene in chrom_coords:
			if ID == gene.internalID:
				transcripts.append(gene)
		return transcripts
	
	def keys(self):
		"""
		Return a list of the keys - duplicating the function of a dictionary
		"""
		return self.gene_coords.keys()

	def getNumGenes(self):
		"""
		Return the number of transcripts
		"""
		num = 0;
		for c in self.gene_coords.keys():
			num += len(self.gene_coords[c]);
		return num;

	def get_strand_specific_genes(self, mystrand):
		"""
		"""
		assert (mystrand == "+" or mystrand =="-")
		mygenes = {}
		for chrom in self.gene_coords.keys():
			mygenes[chrom] = []
			chrom_coords = self.gene_coords[chrom];
			for g in chrom_coords:
				if g.strand == mystrand:
					mygenes[chrom].append(g)
		return mygenes
	
	def output(self, mygenes, outfile):
		o = open(outfile, "w")
		for chrom in mygenes.keys():
			chrom_coords = mygenes[chrom]
			for g in chrom_coords:
				o.write(g.getAll())
		o.close()
		
	def __setitem__(self, name, bedcoord):
		"""
		Sets a new gene coord
		"""
		self.gene_coords[name] = bedcoord


	def __getitem__(self, name):
		"""
		Returns a bed_val indexed by its name or None if no such bed_val exists
		"""
		if self.gene_coords.has_key(name):
			return self.gene_coords[name]
		else: raise UCSCError
    

	def __del__(self):
		"""
		Delete, delete;
		"""
		self.gene_coords.clear()

	def __contains__(self, item):
		"""
		Returns  mapping iterator
		"""
		return self.gene_coords.has_key(item)

	def __iter__(self):
		"""
		Returns mapping iterator
		"""
		return self.gene_coords.iterkeys()

	def __len__(self):
		"""
		Returns number of gene_coords
		"""
		return len(self.gene_coords)

	def __delitem__(self, name):
		"""
		removes a chrom if its name exists in the dictionary
		-- I guess this could possibly be useful at some point
		
		"""
		if self.gene_coords.has_key(name):
			del self.gene_coords[name]
		else: raise UCSCError

def getBoundariesMergedTranscripts(transcripts):
	# need to make sure they are on the same strand from outside
	starts = [transcript.txStart for transcript in transcripts]
	ends = [transcript.txEnd for transcript in transcripts]
	start = min(starts)
	end = max(ends)
	return [(start, end)]

def getMergedExonicRegions(transcripts):
	# Return the merged exons in the format of list of (start, end)
	all_exons = []
	for transcript in transcripts:
		all_exons += transcript.getExons()
	all_exons = sorted(all_exons, key = itemgetter(0))
	return Utility_extended.union(all_exons)


def getSharedExonicRegions(transcripts, min_width=5):
	"""
	shared_exons: a list of (start,end), which are shared among all transcripts in input, might breaking up existing exons.
	sorted
	"""
	
	shared_exons = transcripts[0].getExons()
	for index in range(1, len(transcripts)):
		current_exons = transcripts[index].getExons()
		shared_exons = Utility_extended.intersect(shared_exons, current_exons, min_width)
	return shared_exons #sorted

def getInternalSharedExonicRegions(transcripts, min_width=5):
	"""
	This is used to avoid the interference of promoters and TES.
	"""
	shared_exons = getSharedExonicRegions(transcripts, min_width)
	if len(shared_exons) <= 2:
		return []
	else: 
		return shared_exons[1:-1]
		
def getSharedIntronicRegions(transcripts, min_width=5):
	"""
	shared_introns: a list of (start,end), which are shared among all transcripts in input, might breaking up existing introns.
	sorted
	"""

	shared_introns = transcripts[0].getIntrons()
	for index in range(1, len(transcripts)):
		current_introns = transcripts[index].getIntrons()
		#if len(current_introns) == 0 and len(shared_introns) == 0:
			#for transcript in transcripts:
				#print transcript.getAll()
		shared_introns = Utility_extended.intersect(shared_introns, current_introns, min_width)
	return shared_introns #sorted

def getShared3UTR(transcripts, min_width=5):
	my_3UTRs=[]
	for transcript in transcripts:
		my_3UTR = transcript.get3UTR(0,0)
		if len(my_3UTR) > 0:
			my_3UTRs.append(my_3UTR)
	return Utility_extended.shared(my_3UTRs)

def getShared5UTR(transcripts, min_width=5):
	my_5UTRs=[]
	for transcript in transcripts:
		my_5UTR = transcript.get5UTR(0,0)
		if len(my_5UTR) > 0:
			my_5UTRs.append(my_5UTR)
	return Utility_extended.shared(my_5UTRs)
