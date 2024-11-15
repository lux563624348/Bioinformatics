import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import BED
import GenomeData
import make_graph_file
import SeparateByChrom
import multiprocessing
import get_total_tag_counts


from numpy import *
import scipy
from scipy import *
from scipy import stats
from scipy import special
from scipy import optimize

plus = re.compile("\+");
minus = re.compile("\-");

def get_bed_coordinates(file,chroms,chrom_length,fragment_size):
	infile = open(file);
	postive_tag_counts = 0.0;
	negative_tag_counts = 0.0;
	shift = int(round(fragment_size/2));
	taglist = {};

	for chrom in chroms:
		taglist[chrom]=[]
	
	for line in infile:
		""" check to make sure not a header line """
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			chrom=sline[0]
			if chrom in chroms:
				if atoi(sline[2]) < chrom_length[chrom] and atoi(sline[1]) > 0:
					if plus.match(sline[5]):
						position = atoi(sline[1]) + shift; 
						#If the position is beyond limit then don't shift.
						if position >= chrom_length[chrom]: 
							position = chrom_length[chrom]-1; 
						taglist[chrom].append(position);
						postive_tag_counts += 1.0;
					elif minus.match(sline[5]):
						position = atoi(sline[2]) - 1 - shift;
						# in case the shift move the positions
						# beyond zero, use zero
						if position < 0: position = 0; #UCSC genome coordinate is 0-based
						taglist[chrom].append(position);
						negative_tag_counts += 1.0;
	for chrom in chroms:
		taglist[chrom].sort();

	infile.close();
	return taglist;

def make_histogram(taglist, chrom, chrom_length, window_size):
	bed_vals = {};
	
	if(len(taglist)>0):
		current_window_start = (taglist[0]/window_size)*window_size; 
		tag_count_in_current_window = 1;
		for i in range(1, len(taglist)):
			start = (taglist[i]/window_size)*window_size;		
			if start == current_window_start: tag_count_in_current_window += 1;
			elif start > current_window_start:
				# All the tags in the previous window have been counted 
				current_window_end = current_window_start + window_size -1;
				# if the window goes beyond the chromsome limit, it is discarded. 
				if current_window_end < chrom_length:
					bed_vals[current_window_start]= tag_count_in_current_window;
					# write the window to file
				current_window_start = start;
				tag_count_in_current_window = 1;
			else:
				print 'Something is wrong!!!!!!!';
				
		current_window_end = current_window_start + window_size -1;
		# if the window goes beyond the chromsome limit, it is discarded. 
		if current_window_end < chrom_length:	
			bed_vals[current_window_start]= tag_count_in_current_window;
			
	return bed_vals;
 
	
def calculate_cost(binwidth,taglist,chrom,chrom_length):
	binwidth=int(binwidth)
	if (binwidth < 1):
		return 100000000.0;
	histogram=make_histogram(taglist,chrom,chrom_length,binwidth)
	#print "chrom_length: ", chrom_length, " binwidth: ",binwidth;
	vals=histogram.values();
	cost=(2*mean(vals)-var(vals))/(pow(binwidth,2))
	return cost

if __name__ == "__main__":

	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
					dest="species", help="mm8,hg18,dm2,etc", metavar="<str>")
	parser.add_option("-b", "--bed_file", action="store", type="string",
					dest="bedfile", help="bed file to make graph file of",
					metavar="<file>")
	parser.add_option("-f", "--fragment_size", action="store", type="int",
					dest="fragsize", help="fragment size for chipseq experiment",
					metavar="<file>")
	(opt, args) = parser.parse_args(sys.argv)

	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species];
	genome_length=float(sum(chrom_lengths.values()))
	total_read_count = get_total_tag_counts.get_total_tag_counts(opt.bedfile);
	tag_coords=get_bed_coordinates(opt.bedfile,chroms,chrom_lengths,opt.fragsize)
			
	#print total_read_count
	#print "Genome size: " , genome_length*0.76

	pool=multiprocessing.Pool()
	results={}
	binwidths=[]

	for chrom in chroms:
		results[chrom]=pool.apply_async(optimize.fminbound,(calculate_cost,50,1.e6,(tag_coords[chrom],chrom,chrom_lengths[chrom],),))
		binwidth=results[chrom].get()
		if (int(binwidth) > 0):
			binwidths.append(binwidth)

	pool.close()
	pool.join()

	window_size=int(round(median(binwidths),-1))
#		    result=optimize.fminbound(calculate_cost,50,1.e6,args=(tag_coords[chrom],chrom,chrom_lengths[chrom],))
#		    break
	print "Windowsize: ",window_size
	# count tags, read starts per chromosome
	# pass tags into function per chromosome. Function calculates cost(windowsize).
	# in function, create graphs depending on windowsize.
	# calculate mean and variance for graphs
	# calculate Cost.
	# minimize function
	    
	  
