import os
import sys
import numpy as np
import pandas as pd
import collections
import itertools
from mirnylib import h5dict, genome
from hiclib import mapping, fragmentHiC
import cooler
from multiprocessing import Pool

genome_db = genome.Genome('/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa')

# this part has run already
'''
fname_list = ['DKO_CD8_p%d_mapped_reads.hdf5' % i for i in range(1,7)]

mapped_reads_list = []
for fname in fname_list:
    mapped_reads_list.append(h5dict.h5dict(fname))

combine_mapped_reads = h5dict.h5dict('DKO_CD8_mapped_reads.hdf5')

combine_mapped_reads['chrms1'] = np.hstack(map(lambda k: k['chrms1'], mapped_reads_list))
combine_mapped_reads['chrms2'] = np.hstack(map(lambda k: k['chrms2'], mapped_reads_list))
combine_mapped_reads['cuts1'] = np.hstack(map(lambda k: k['cuts1'], mapped_reads_list))
combine_mapped_reads['cuts2'] = np.hstack(map(lambda k: k['cuts2'], mapped_reads_list))
combine_mapped_reads['misc'] = map(lambda k: k['misc'], mapped_reads_list)[0]
combine_mapped_reads['strands1'] = np.hstack(map(lambda k: k['strands1'], mapped_reads_list))
combine_mapped_reads['strands2'] = np.hstack(map(lambda k: k['strands1'], mapped_reads_list))

fragments = fragmentHiC.HiCdataset(
    filename='DKO_CD8_fragment_dataset.hdf5',
    enzymeName='HindIII',
    genome=genome_db,
    mode='w')

fragments.parseInputData(
    dictLike='DKO_CD8_mapped_reads.hdf5', keepSingleSided=False, keepSameFragment=True, keepReadsMolecules=True)

fragments.filterDuplicates()

fragments._sortData()

'''
# run ICE

fragments = fragmentHiC.HiCdataset('DKO_CD8_fragment_dataset.hdf5', genome_db, enzymeName='HindIII')
fragments.saveCooler('DKO_CD8-res-10K.cool', 10000)

filename = 'DKO_CD8-res-10K.cool'
mypool = Pool(8)

mycooler = cooler.Cooler(filename)
N = mycooler.info["nnz"]
Nbin = mycooler.info["nbins"]
if N / Nbin < 30:
    print("couldn't correct {0}".format(filename))
    print("N/Nbin={}".format(N/Nbin))
    raise 

print("Starting ICE", filename); sys.stdout.flush()

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
    chunk.loc[:,['chrom1', 'mid1', 'chrom2', 'mid2', 'count']].to_csv('fit-hi-c/DKO_CD8_contactCounts.txt', index=None, header=None, sep='\t', mode='a')

bin_df = mycooler.bins()[:]
bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)
bin_df['extraField'] = 0
bin_df['marginalContactCount'] = bin_df.apply(lambda row: 0 if np.isnan(row['weight']) else 1, axis=1)
bin_df['mappable'] = 0
bin_df.loc[:,['chrom', 'extraField', 'mid', 'marginalContactCount', 'mappable']].to_csv('fit-hi-c/DKO_CD8_marginal_contactCounts.txt', index=None, header=None, sep='\t')

bin_df['bias'] =  bin_df['weight'] / bin_df['weight'].mean()
bin_df.loc[np.isnan(bin_df.bias), 'bias'] = 0
bin_df.loc[:,['chrom', 'mid', 'bias']].to_csv('fit-hi-c/DKO_CD8_ICE_bias.txt', index=None, header=None, sep='\t')

# generate Juicebox input

for chr1, chr2 in itertools.product(mycooler.chromnames, mycooler.chromnames):
    chunk = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chr1, chr2).dropna()
    chunk['mid1'] = ((chunk['start1'] + chunk['end1']) / 2).astype(int)
    chunk['mid2'] = ((chunk['start2'] + chunk['end2']) / 2).astype(int)
    chunk.loc[:,['chrom1', 'mid1', 'chrom2', 'mid2', 'balanced']].to_csv('Juicebox/DKO_CD8_Juicebox_input.txt', sep='\t', index=None, header=None, mode='a')
    chunk['str1'] = 0
    chunk['frag1'] = 0
    chunk['str2'] = 0
    chunk['frag2'] = 1
    chunk.loc[:,['str1', 'chrom1', 'mid1', 'frag1', 'str2', 'chrom2', 'mid2', 'frag2', 'balanced']].to_csv('Juicebox/DKO_CD8_Juicebox_input.txt', sep='\t', index=None, header=None, mode='a')
