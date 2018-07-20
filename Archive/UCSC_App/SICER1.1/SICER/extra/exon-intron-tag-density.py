#!/usr/bin/env python
# Copyright (c) 2009 GWU & NHLBI, NIH
# Authors: Chongzhi Zang, Dustin Schones, Weiqun Peng and Keji Zhao
#
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
# zang@gwmail.gwu.edu).


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import bisect

import BED
import UCSC
import GenomeData
import Utility
import SeparateByChrom
import associate_tags_with_regions

plus = re.compile("\+");
minus = re.compile("\-");



def getExonIntronDensities(gene_coords, bed_vals, num_steps, fragment_size):
	shift = int(round(fragment_size/2))
	exon_sums = [];
	intron_sums = [];
	exon_sizes = [];
	intron_sizes = [];
	for i in range(int(num_steps)):
		exon_sums.append(0);
		intron_sums.append(0);
		exon_sizes.append(0);
		intron_sizes.append(0);
	for chrom in bed_vals.keys():
		tag_positions = [];
		for item in bed_vals[chrom]:
			if plus.match(item.strand):
				tag_positions.append(item.start + shift);
			elif minus.match(item.strand):
				tag_positions.append(item.end - 1 - shift);
		tag_positions.sort();
            
		if chrom in gene_coords.keys() and len(tag_positions)>0:
			for g in gene_coords[chrom]:
				#print g.name
				steps = min(num_steps, g.exonCount-1)
				exon_Starts = g.exonStarts.split(',')
				exon_Ends = g.exonEnds.split(',')
				assert len(exon_Starts) == len(exon_Ends)
				if steps > 0:
					if plus.match(g.strand):
						# exon 0, 1, 2, ...g.exonCount-2 or num_steps-1
						for i in range(0, int(steps)):
							exon_sums[i] += associate_tags_with_regions.countTagsInWindow(int(exon_Starts[i]), int(exon_Ends[i]), tag_positions);
							intron_sums[i] += associate_tags_with_regions.countTagsInWindow(int(exon_Ends[i]), int(exon_Starts[i+1]), tag_positions);
							exon_sizes[i] += abs(int(exon_Ends[i]) - int(exon_Starts[i]));
							intron_sizes[i] += abs(int(exon_Starts[i+1]) - int(exon_Ends[i]));
						# last exon
						if ((g.exonCount-1) < int(num_steps)):
							exon_sums[g.exonCount-1] += associate_tags_with_regions.countTagsInWindow(int(exon_Starts[g.exonCount-1]), int(exon_Ends[g.exonCount-1]), tag_positions);
							exon_sizes[g.exonCount-1] += abs(int(exon_Ends[g.exonCount-1]) - int(exon_Starts[g.exonCount-1]));
					elif minus.match(g.strand):
						for i in range(0, int(steps)):
						#print exon_Starts
							exon_sums[i] += associate_tags_with_regions.countTagsInWindow(int(exon_Starts[-2-i]), int(exon_Ends[-2-i]), tag_positions);
							intron_sums[i] += associate_tags_with_regions.countTagsInWindow(int(exon_Ends[-3-i]), int(exon_Starts[-2-i]), tag_positions);
							exon_sizes[i] += abs(int(exon_Ends[-2-i]) - int(exon_Starts[-2-i]));
							intron_sizes[i] += abs(int(exon_Starts[-2-i]) - int(exon_Ends[-3-i]));
						if ((g.exonCount-1) < int(num_steps)):
							exon_sums[g.exonCount-1] += associate_tags_with_regions.countTagsInWindow(int(exon_Starts[0]), int(exon_Ends[0]), tag_positions);
							exon_sizes[g.exonCount-1] += abs(int(exon_Ends[0]) - int(exon_Starts[0]));
	return exon_sums, intron_sums, exon_sizes, intron_sizes;




def main(argv):
	parser = OptionParser()
	parser.add_option("-k", "--known_gene_file", action="store", type="string",
			dest="genefile", help="file with known gene info", metavar="<file>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bedfile", help="file with tags in bed format", metavar="<file>")
	parser.add_option("-n", "--name", action="store", type="string",
			dest="name", help="name for plotting", metavar="<str>")
	parser.add_option("-s", "--species", action="store", type="string",
			dest="species", help="species", metavar="<str>")
	parser.add_option("-c", "--number", action="store", type="int",
			dest="exonnumber", help="number of exons", metavar="<int>")
	parser.add_option("-f", "--fragmentsize", action="store", type="int",
			dest="fragment_size", help="fragment size", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		chrom_length = GenomeData.species_chrom_lengths[opt.species];
		
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	

    	gene_coords = UCSC.KnownGenes(opt.genefile);
	
	#separate by chrom
	if Utility.fileExists(opt.bedfile):
		SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1');
	else:
		print opt.beddfile, " not found";
		sys.exit(1)
	
	total_exon_sums = [0]*opt.exonnumber;
	total_intron_sums = [0]*opt.exonnumber;
	total_exon_sizes = [0]*opt.exonnumber;
	total_intron_sizes = [0]*opt.exonnumber;
	total_num_tags = 0;
	
	for chrom in chroms:
		read_file = chrom + ".bed1";
		bed_vals = BED.BED(opt.species, read_file, "BED6", 0);
    		total_num_tags += bed_vals.getNumVals();
		(exon_counts, intron_counts, exon_seq_sizes, intron_seq_sizes) = getExonIntronDensities(gene_coords, bed_vals, opt.exonnumber, opt.fragment_size);
  		for j in range(opt.exonnumber):
			total_exon_sums[j] += exon_counts[j];
			total_intron_sums[j] += intron_counts[j];
			total_exon_sizes[j] += exon_seq_sizes[j];
			total_intron_sizes[j] += intron_seq_sizes[j];

	""" print everything out to a file """
    	outfilename = '%s-exon-intron-scores' % opt.name;
    	outFile = open(outfilename, 'w');    	
	for j in range(len(exon_counts)):
		exon_density = float(total_exon_sums[j]) / float(total_exon_sizes[j]);
		intron_density = float(total_intron_sums[j]) / float(total_intron_sizes[j]);
		exon_density /= float(total_num_tags);
		intron_density  /= float(total_num_tags);
		outline = str(j+1) + " " + str(exon_density) + " " + str(intron_density) + "\n";
		outFile.write(outline);
	outFile.close();

	SeparateByChrom.cleanup(chroms,'.bed1');

if __name__ == "__main__":
	main(sys.argv)


        
