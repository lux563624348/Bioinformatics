import os
import sys
import numpy as np
import pandas as pd
import collections
from mirnylib import h5dict, genome, plotting
from hiclib import mapping, fragmentHiC, highResBinnedData

genome_db = genome.Genome('/home/xli/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa')

BD = highResBinnedData.HiResHiC(genome_db, 10000, 'WT_CD8_heatmap-res-10K-IC.hdf5')
BD.loadData('WT_CD8_heatmap-res-10K.hdf5')


BD.removeDiagonal()

BD.iterativeCorrection()



outdir = '/home/xli/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Combine_Jun_Jul_2016/iterative_mapping_3/3d_genome_browser/3d_genome_browser_WT_CD8'

pd.DataFrame({'chrmLabels': BD.genome.chrmLabels, 'chrmLens': BD.genome.chrmLens}).to_csv(outdir + 'mm9_genome_size.txt', index=None, header=None, sep='\t')
pd.DataFrame({'chrmLabels': BD.genome.chrmLabels, 'filename': map(lambda chrm: outdir + 'WT_CD8_' + chrm + '.matrix', BD.genome.chrmLabels)}).to_csv(outdir + 'WT_CD8_intrachrom_matrix_list.txt', index=None, header=None, sep='\t')

for index, chrmLabel in BD.genome.idx2label.iteritems():
    key = '(%d, %d)' % (index, index)
    outputfile = 'WT_CD8_' +  chrmLabel + '.matrix'
    np.savetxt(outputfile, BD._h5dict[key], fmt='%.2e', delimiter='\t')
