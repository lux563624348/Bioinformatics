#!/usr/bin/env python
# Authors: Dustin E Schones, Chongzhi Zang, Weiqun Peng and Keji Zhao
# Disclaimer
# 
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import GenomeData

import HTSeq
import itertools

def get_read_length(alt_iv_seq):
        read_length = 0
        for alt_iv in alt_iv_seq:
                read_length += alt_iv.length
        return read_length

def reverse_strand( iv ):
        if iv.strand == "+":
                iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
        elif iv.strand == "-":
                iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
        return iv_reversed

def combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq):
        combine_alt_iv_seq = []
        alt_iv_seq = alt_first_iv_seq + alt_second_iv_seq
        alt_iv_seq.sort(key=lambda x: x.start)
        combine_alt_iv_seq.append(alt_iv_seq[0].copy())
        for alt in itertools.islice(alt_iv_seq, 1, None):
                if alt.overlaps(combine_alt_iv_seq[-1]):
                        combine_alt_iv_seq[-1].extend_to_include(alt)
                else:
                        combine_alt_iv_seq.append(alt.copy())
        return combine_alt_iv_seq

def makeGraphFile(chroms, chrom_lengths, window_size, bamfile, outfile):
        ga = HTSeq.GenomicArray(chroms, stranded=False, typecode='d')
        
        bam_reader = HTSeq.BAM_Reader(bamfile)
        for alt_first, alt_second in HTSeq.pair_SAM_alignments(bam_reader):
                if alt_first == None or alt_second == None:
                        continue
                if alt_first.aligned and alt_first.optional_field("NH") == 1 and alt_second.aligned and alt_second.optional_field("NH") == 1:
                        if alt_first.iv.chrom != alt_second.iv.chrom or alt_first.iv.strand == alt_second.iv.strand or alt_first.iv.chrom not in chroms:
                                continue
                        
                        alt_first_iv_seq = [ co.ref_iv for co in alt_first.cigar if co.type == "M" and co.size > 0 ]
                        alt_second_iv_seq = [ reverse_strand(co.ref_iv) for co in alt_second.cigar if co.type == "M" and co.size > 0 ]                            
                        alt_iv_seq = combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)
                        
                        read_length = get_read_length(alt_iv_seq)
                        for alt_iv in alt_iv_seq:
                                ga[alt_iv] += 1.0 / read_length
        
        with open(outfile, 'w') as f:
                for chrom in chroms:
                        chrom_length = chrom_lengths[chrom]
                        num_windows = chrom_length / window_size                       
                        for i in range(num_windows):
                                count_in_window = 0
                                window_start = i * window_size
                                window_end = (i+1) * window_size
                                window_iv = HTSeq.GenomicInterval(chrom, window_start, window_end)
                                for iv, value in ga[window_iv].steps():
                                        count_in_window += value * iv.length
                                count_in_window = int(count_in_window)
                                if count_in_window != 0:
                                        outline = chrom + '\t' + str(window_start) + '\t' + str(window_end) + '\t' + str(count_in_window) + '\n'
                                        f.write(outline)	

def makeIslandFilteredGraphFile(chroms, chrom_lengths, window_size, bamfile, islandbedfile, outfile):
        ga = HTSeq.GenomicArray(chroms, stranded=False, typecode='d')
        
        bam_reader = HTSeq.BAM_Reader(bamfile)
        for alt_first, alt_second in HTSeq.pair_SAM_alignments(bam_reader):
                if alt_first == None or alt_second == None:
                        continue
                if alt_first.aligned and alt_first.optional_field("NH") == 1 and alt_second.aligned and alt_second.optional_field("NH") == 1:
                        if alt_first.iv.chrom != alt_second.iv.chrom or alt_first.iv.strand == alt_second.iv.strand or alt_first.iv.chrom not in chroms:
                                continue
                        
                        alt_first_iv_seq = [ co.ref_iv for co in alt_first.cigar if co.type == "M" and co.size > 0 ]
                        alt_second_iv_seq = [ reverse_strand(co.ref_iv) for co in alt_second.cigar if co.type == "M" and co.size > 0 ]                            
                        alt_iv_seq = combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)
                        
                        read_length = get_read_length(alt_iv_seq)
                        for alt_iv in alt_iv_seq:
                                ga[alt_iv] += 1.0 / read_length
                                
        ga_island = HTSeq.GenomicArray(chroms, stranded=False, typecode='d')
        bedfile = HTSeq.BED_Reader(islandbedfile)
        for alt in bedfile:
                for iv, value in ga[alt.iv].steps():
                        ga_island[iv] += value
        
        with open(outfile, 'w') as f:
                for chrom in chroms:
                        chrom_length = chrom_lengths[chrom]
                        num_windows = chrom_length / window_size                       
                        for i in range(num_windows):
                                count_in_window = 0
                                window_start = i * window_size
                                window_end = (i+1) * window_size
                                window_iv = HTSeq.GenomicInterval(chrom, window_start, window_end)
                                for iv, value in ga_island[window_iv].steps():
                                        count_in_window += value * iv.length
                                count_in_window = int(count_in_window)
                                if count_in_window != 0:
                                        outline = chrom + '\t' + str(window_start) + '\t' + str(window_end) + '\t' + str(count_in_window) + '\n'
                                        f.write(outline)        

def main(argv):
        """
        Note the window_size and the fragment_size are both input as strings, as they are used in
        a shell script in makeGraphFile.
        """
        parser = OptionParser()
        parser.add_option("-s", "--species", action="store", type="string",
                          dest="species", help="mm8,hg18,dm2,etc", metavar="<str>")
        parser.add_option("-b", "--bam_file", action="store", type="string",
                          dest="bamfile", help="bam file to make graph file of",
                          metavar="<file>")
        parser.add_option("-w", "--window_size", action="store", type="int",
                          dest="window_size", help="window size", metavar="<int>")
        parser.add_option("-i", "--islandfile", action="store", type="string",
                          dest="islandbedfile", metavar="<file>",
                          help="island file")        
        parser.add_option("-o", "--outfile", action="store", type="string",
                          dest="outfile", help="output bed summary file name",
                          metavar="<file>")
        (opt, args) = parser.parse_args(argv)
        if len(argv) < 8:
                parser.print_help()
                sys.exit(1)

        if opt.species in GenomeData.species_chroms.keys():
                chroms = GenomeData.species_chroms[opt.species];
                chrom_lengths = GenomeData.species_chrom_lengths[opt.species];      
                if opt.islandbedfile:
                        makeIslandFilteredGraphFile(chroms, chrom_lengths, opt.window_size, opt.bamfile, opt.islandbedfile, opt.outfile)
                else:
                        makeGraphFile(chroms, chrom_lengths, opt.window_size, opt.bamfile, opt.outfile)
        else:
                print opt.species + " is not in the species list ";



if __name__ == "__main__":
        main(sys.argv)



