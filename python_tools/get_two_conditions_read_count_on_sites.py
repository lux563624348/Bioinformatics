import HTSeq
import sys, os, re
from string import *
from optparse import OptionParser
import numpy
import collections
import itertools

def main(argv):
    parser = OptionParser()
    parser.add_option("-a", "--tags_bed_file_C1", action="store", type="string", dest="tagsfile_C1", metavar="<file>", help="input ChIP-seq tags bed file in condition 1")
    parser.add_option("-b", "--tags_bed_file_C2", action="store", type="string", dest="tagsfile_C2", metavar="<file>", help="input ChIP-seq tags bed file in condition 2")
    parser.add_option("-p", "--peaks_file", action="store", type="string", dest="peaksfile", metavar="<file>", help="peaks file")	
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="fragment size determines the shift (half of fragment_size of ChIP-seq read position, in bps")
    parser.add_option("-s", "--output_summary_file", action="store", type="string", dest="summaryfile", metavar="<file>", help="output summary file for read count in two conditions")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    fragment_size = opt.fragment_size

    chip_C1 = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    tagsfile_C1 = HTSeq.BED_Reader(opt.tagsfile_C1)
    for alt in tagsfile_C1:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        chip_C1[alt_pos] += 1

    chip_C2 = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    tagsfile_C2 = HTSeq.BED_Reader(opt.tagsfile_C2)
    for alt in tagsfile_C2:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        chip_C2[alt_pos] += 1	

    summaryfile = open(opt.summaryfile, "w")
    header = "chrom" + "\t" + "start" + "\t" + "end" + "\t" + "tags for C1" + "\t" + "tags for C2" + "\n"
    summaryfile.write(header)

    peaksfile = HTSeq.BED_Reader(opt.peaksfile)
    for peak in peaksfile:
        peaks_iv_count_C1 = 0
        for count in chip_C1[peak.iv]:
            peaks_iv_count_C1 += count
        peaks_iv_count_C2 = 0
        for count in chip_C2[peak.iv]:
            peaks_iv_count_C2 += count  
        
        summaryline = peak.iv.chrom + "\t" + str(peak.iv.start) + "\t" + str(peak.iv.end) + "\t" + str(peaks_iv_count_C1) + "\t" + str(peaks_iv_count_C2) + "\n"
        summaryfile.write(summaryline)        
        
    summaryfile.close()


if __name__ == "__main__":
    main(sys.argv)


