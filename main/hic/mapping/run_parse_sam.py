import sys, os
import numpy as np
from mirnylib import h5dict, genome
from hiclib_old import mapping, fragmentHiC, binnedData
from optparse import OptionParser

def main(argv):
    parser = OptionParser()	
    parser.add_option("-g", "--genome_path", action="store", type="string", dest="genome_path", metavar="<file>", help="path for genome.fa")
    parser.add_option("-1", "--bamfile_R1", action="store", type="string", dest="bamfile_R1", metavar="<file>", help="R1 bam file")
    parser.add_option("-2", "--bamfile_R2", action="store", type="string", dest="bamfile_R2", metavar="<file>", help="R2 bam file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output hdf5 file")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)


    genome_db = genome.Genome(opt.genome_path)

    mapped_reads = h5dict.h5dict(opt.outfile)

    mapping.parse_sam(
        sam_basename1=opt.bamfile_R1,
        sam_basename2=opt.bamfile_R2,
        out_dict=mapped_reads,
        genome_db=genome_db)
    
if __name__ == "__main__":
    main(sys.argv)
