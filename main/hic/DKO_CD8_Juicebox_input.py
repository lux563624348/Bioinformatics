import os
import sys
import numpy as np
import pandas as pd
from mirnylib import h5dict, genome, plotting
from hiclib import mapping, fragmentHiC, highResBinnedData

genome_db = genome.Genome('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa')

fragments = fragmentHiC.HiCdataset('DKO_CD8_fragment_dataset.hdf5', genome_db, enzymeName='HindIII')

chrms1_array = np.empty_like(fragments.h5dict['chrms1'], dtype='|S5')
for index, chrmLabel in fragments.genome.idx2label.iteritems():
    chrms1_array[np.where(fragments.h5dict['chrms1'] == index)] = chrmLabel

chrms2_array = np.empty_like(fragments.h5dict['chrms2'], dtype='|S5')
for index, chrmLabel in fragments.genome.idx2label.iteritems():
    chrms2_array[np.where(fragments.h5dict['chrms2'] == index)] = chrmLabel

k = 10
index_list = np.linspace(0, fragments.N, k + 1, dtype=int)

for i in range(1,k+1):
    print i; sys.stdout.flush()
    start, end = index_list[i-1], index_list[i]
    df = pd.DataFrame({'str1': 1 - fragments.h5dict['strands1'][start : end],
                       'chr1': chrms1_array[start : end],
                       'pos1': fragments.h5dict['cuts1'][start : end],
                       'frag1': 0,
                       'str2': 1 - fragments.h5dict['strands2'][start : end],
                       'chr2': chrms2_array[start : end],
                       'pos2': fragments.h5dict['cuts2'][start : end],
                       'frag2': 1
                      },
                      columns=['str1','chr1','pos1','frag1','str2','chr2','pos2','frag2'])

    df.to_csv('Juicebox/DKO_CD8_Juicebox_input.txt', sep='\t', index=None, header=None, mode='a')
