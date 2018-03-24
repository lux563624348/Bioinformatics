import sys, os
import re
import numpy as np
import pandas as pd
import collections
import itertools
from mirnylib import h5dict, genome
from hiclib import mapping, fragmentHiC
import cooler
from multiprocessing import Pool
from optparse import OptionParser
import subprocess

def main(argv):
    parser = OptionParser()	
    parser.add_option("-g", "--genome_path", action="store", type="string", dest="genome_path", metavar="<file>", help="path for genome.fa")
    parser.add_option("-n", "--name", action="store", type="string", dest="name", metavar="<file>", help="sample name")
    parser.add_option("-o", "--outdir", action="store", type="string", dest="outdir", metavar="<file>", help="output directory")
    parser.add_option("-r", "--resolution", action="store", type="int", dest="resolution", metavar="<int>", help="heatmap resolution")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)
        
        
    # output name
    reads_hdf5_file = os.path.join(opt.outdir, '%s_mapped_reads.hdf5' % opt.name)
    fragment_hdf5_file = os.path.join(opt.outdir, '%s_fragment_dataset.hdf5' % opt.name)
    cool_file = os.path.join(opt.outdir, '%s-res-%dK.cool' % (opt.name, opt.resolution / 1000))
    
    if not os.path.exists(os.path.join(opt.outdir, 'fit-hi-c/')):
        os.makedirs(os.path.join(opt.outdir, 'fit-hi-c/'))
        
    fit_hic_contactCounts_file = os.path.join(opt.outdir, 'fit-hi-c/%s_contactCounts.txt' % opt.name)
    fit_hic_marginal_contactCounts_file = os.path.join(opt.outdir, 'fit-hi-c/%s_marginal_contactCounts.txt' % opt.name)
    fit_hic_ICE_bias_file = os.path.join(opt.outdir, 'fit-hi-c/%s_ICE_bias.txt' % opt.name)
    
    if not os.path.exists(os.path.join(opt.outdir, 'Juicebox/')):
        os.makedirs(os.path.join(opt.outdir, 'Juicebox/'))    
    
    juicebox_file = os.path.join(opt.outdir, 'Juicebox/%s_Juicebox_input.txt' % opt.name)
      

    genome_db = genome.Genome(opt.genome_path)
    
    # run ICE
    
    fragments = fragmentHiC.HiCdataset(fragment_hdf5_file, genome_db, enzymeName='HindIII')
    fragments.saveCooler(cool_file, opt.resolution)
    
    mypool = Pool(8)
    
    mycooler = cooler.Cooler(cool_file)
    N = mycooler.info["nnz"]
    Nbin = mycooler.info["nbins"]
    if N / Nbin < 30:
        print("couldn't correct {0}".format(cool_file))
        print("N/Nbin={}".format(N/Nbin))
        raise 
    
    print("Starting ICE", cool_file); sys.stdout.flush()
    
    with mycooler.open() as h5:
        
        bias, stats = cooler.ice.iterative_correction(h5, chunksize= min(N//8+5, 40000000), map=mypool.map, 
                                                      min_nnz=10, use_lock=False, mad_max=4, ignore_diags=2, rescale_marginals=False,
                                                      max_iters=500)
        
    with mycooler.open('r+') as h5:
        if 'weight' in h5['bins']:
            del h5["bins/weight"]
    
        h5['bins'].create_dataset('weight', data=bias)
        h5['bins']['weight'].attrs.update(stats)
        
    mypool.close()
    
    print "ICE done"; sys.stdout.flush()
    
    # generate fit-hi-c input
    
    for chrom in mycooler.chromnames:
        chunk = mycooler.pixels(join=True).fetch(chrom)
        chunk = chunk[(chunk.chrom1 == chunk.chrom2) & (chunk.start1 != chunk.start2)]
        chunk['mid1'] = ((chunk['start1'] + chunk['end1']) / 2).astype(int)
        chunk['mid2'] = ((chunk['start2'] + chunk['end2']) / 2).astype(int)
        chunk.loc[:,['chrom1', 'mid1', 'chrom2', 'mid2', 'count']].to_csv(fit_hic_contactCounts_file, index=None, header=None, sep='\t', mode='a')
    
    bin_df = mycooler.bins()[:]
    bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)
    bin_df['extraField'] = 0
    bin_df['marginalContactCount'] = bin_df.apply(lambda row: 0 if np.isnan(row['weight']) else 1, axis=1)
    bin_df['mappable'] = 0
    bin_df.loc[:,['chrom', 'extraField', 'mid', 'marginalContactCount', 'mappable']].to_csv(fit_hic_marginal_contactCounts_file, index=None, header=None, sep='\t')
    
    bin_df['bias'] =  bin_df['weight'] / bin_df['weight'].mean()
    bin_df.loc[np.isnan(bin_df.bias), 'bias'] = 0
    bin_df.loc[:,['chrom', 'mid', 'bias']].to_csv(fit_hic_ICE_bias_file, index=None, header=None, sep='\t')
    
    # generate Juicebox input
    
    for chr1, chr2 in itertools.product(mycooler.chromnames, mycooler.chromnames):
        chunk = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chr1, chr2).dropna()
        chunk['mid1'] = ((chunk['start1'] + chunk['end1']) / 2).astype(int)
        chunk['mid2'] = ((chunk['start2'] + chunk['end2']) / 2).astype(int)
        chunk['str1'] = 0
        chunk['frag1'] = 0
        chunk['str2'] = 0
        chunk['frag2'] = 1
        chunk.loc[:,['str1', 'chrom1', 'mid1', 'frag1', 'str2', 'chrom2', 'mid2', 'frag2', 'balanced']].to_csv(juicebox_file, sep='\t', index=None, header=None, mode='a')
        
        
    for f in [fit_hic_contactCounts_file, fit_hic_marginal_contactCounts_file, fit_hic_ICE_bias_file, juicebox_file]:
        p = subprocess.Popen('gzip %s' % f, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
        p.wait()
        
        
if __name__ == "__main__":
    main(sys.argv)

    


    