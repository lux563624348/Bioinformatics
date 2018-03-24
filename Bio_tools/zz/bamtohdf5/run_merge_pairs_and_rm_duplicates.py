import sys, os
import re
import numpy as np
from mirnylib import h5dict, genome
from hiclib import mapping, fragmentHiC
from optparse import OptionParser

def main(argv):
    parser = OptionParser()	
    parser.add_option("-g", "--genome_path", action="store", type="string", dest="genome_path", metavar="<file>", help="path for genome.fa")
    parser.add_option("-i", "--indir", action="store", type="string", dest="indir", metavar="<file>", help="input partioned hdf5 files directory")
    parser.add_option("-n", "--name", action="store", type="string", dest="name", metavar="<file>", help="sample name")
    parser.add_option("-o", "--outdir", action="store", type="string", dest="outdir", metavar="<file>", help="output directory")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)
        
    reads_hdf5_file = os.path.join(opt.outdir, '%s_mapped_reads.hdf5' % opt.name)
    fragment_hdf5_file = os.path.join(opt.outdir, '%s_fragment_dataset.hdf5' % opt.name)

    genome_db = genome.Genome(opt.genome_path)

    #fname_list = filter(lambda f: re.search('%s_p\d+_mapped_reads.hdf5' % opt.name, f), os.listdir(opt.indir))

#### Above line is for the filename with p1 p2 p3 case. 

    fname_list = filter(lambda f: re.search('%s_mapped_reads.hdf5' % opt.name, f), os.listdir(opt.indir))



    mapped_reads_list = []
    for fname in fname_list:
        mapped_reads_list.append(h5dict.h5dict(os.path.join(opt.indir, fname)))
    
    combine_mapped_reads = h5dict.h5dict(reads_hdf5_file)
    
    combine_mapped_reads['chrms1'] = np.hstack(map(lambda k: k['chrms1'], mapped_reads_list))
    combine_mapped_reads['chrms2'] = np.hstack(map(lambda k: k['chrms2'], mapped_reads_list))
    combine_mapped_reads['cuts1'] = np.hstack(map(lambda k: k['cuts1'], mapped_reads_list))
    combine_mapped_reads['cuts2'] = np.hstack(map(lambda k: k['cuts2'], mapped_reads_list))
    combine_mapped_reads['misc'] = map(lambda k: k['misc'], mapped_reads_list)[0]
    combine_mapped_reads['strands1'] = np.hstack(map(lambda k: k['strands1'], mapped_reads_list))
    combine_mapped_reads['strands2'] = np.hstack(map(lambda k: k['strands1'], mapped_reads_list))
    
    fragments = fragmentHiC.HiCdataset(
        filename=fragment_hdf5_file,
        enzymeName='HindIII',
        genome=genome_db,
        mode='w')

    fragments.parseInputData(
        dictLike=reads_hdf5_file, keepSingleSided=False, keepSameFragment=True, keepReadsMolecules=True)
    
    fragments.filterDuplicates()
    
    fragments._sortData()

if __name__ == "__main__":
    main(sys.argv)
