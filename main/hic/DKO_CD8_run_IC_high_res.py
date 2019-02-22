import os
import sys
import numpy as np
import pandas as pd
import collections
from mirnylib import h5dict, genome, plotting
from hiclib import mapping, fragmentHiC, highResBinnedData

genome_db = genome.Genome('~/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa')

# this part has run already

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

fragments = fragmentHiC.HiCdataset('DKO_CD8_fragment_dataset.hdf5', genome_db, enzymeName='HindIII')

fragments.saveByChromosomeHeatmap('DKO_CD8_heatmap-res-10K.hdf5', resolution=10000)
fragments.printMetadata(saveTo="DKO_CD8_statistics.txt")
#'''
#

BD = highResBinnedData.HiResHiC(genome_db, 10000, 'DKO_CD8_heatmap-res-10K-IC.hdf5')
BD.loadData('DKO_CD8_heatmap-res-10K.hdf5')

BD.removeDiagonal()

bias_chrms = []
bias_fragmentMids = []

mid_point_dict = collections.defaultdict(lambda: {})
resolution = BD.resolution
for i in range(BD.chromosomeCount):
    last_index = BD.genome.chrmLensBin[i]
    for j in range(BD.genome.chrmLensBin[i]):
        if j != last_index:
            mid_point_dict[i][j] = j * resolution + resolution / 2
        else:
            mid_point_dict[i][j] = (j * resolution + BD.genome.chrmLens[i]) / 2

        bias_chrms.append(BD.genome.idx2label[i])
        bias_fragmentMids.append(mid_point_dict[i][j])

for index, chrmLabel in BD.genome.idx2label.iteritems():
    marginal_array = BD.getMarginals()[index]
    num_entries = len(marginal_array)
    f1_df = pd.DataFrame({'chr': [chrmLabel] * num_entries,
                          'extraField': 0,
                          'fragmentMid': map(lambda i: mid_point_dict[index][i], range(num_entries)),
                          'contactCount': marginal_array,
                          'mappable': 0
                          },
			  columns=['chr','extraField','fragmentMid','contactCount','mappable'])
    f1_df['contactCount'] = f1_df['contactCount'].astype(int)
    f1_df.to_csv('fit-hi-c/DKO_CD8_marginal_contactCounts.txt', index=None, header=None, sep='\t', mode='a')

for index, chrmLabel in BD.genome.idx2label.iteritems():
    print chrmLabel; sys.stdout.flush()
    key = '(%d, %d)' % (index, index)
    matrix = np.triu(BD._h5dict[key])
    num_entries = len(np.nonzero(matrix)[0])
    f2_df = pd.DataFrame({'chr1': [chrmLabel] * num_entries,
                          'fragmentMid1': map(lambda i: mid_point_dict[index][i], np.nonzero(matrix)[0]),
                          'chr2': [chrmLabel] * num_entries,
                          'fragmentMid2': map(lambda i: mid_point_dict[index][i], np.nonzero(matrix)[1]),
                          'contactCount': matrix[np.nonzero(matrix)]
                         },
                         columns=['chr1','fragmentMid1','chr2','fragmentMid2','contactCount'])
    f2_df['contactCount'] = f2_df['contactCount'].astype(int)
    f2_df.to_csv('fit-hi-c/DKO_CD8_contactCounts.txt', sep='\t', index=None, header=None, mode='a')


BD.iterativeCorrection()

bias_df = pd.DataFrame({'chr': bias_chrms,
                        'fragmentMid': bias_fragmentMids,
                        'bias': BD.biases},
                        columns = ['chr', 'fragmentMid', 'bias'])

bias_df.to_csv('fit-hi-c/DKO_CD8_ICE_bias.txt', index=None, header=None, sep='\t')

for index, chrmLabel in BD.genome.idx2label.iteritems():
    print chrmLabel; sys.stdout.flush()
    key = '(%d, %d)' % (index, index)
    matrix = np.triu(BD._h5dict[key])
    num_entries = len(np.nonzero(matrix)[0])
    Juicebox_df = pd.DataFrame({'str1': [0] * num_entries,
                                'chr1': [chrmLabel] * num_entries,
                                'fragmentMid1': map(lambda i: mid_point_dict[index][i], np.nonzero(matrix)[0]),
                                'frag1': [0] * num_entries,
                                'str2': [0] * num_entries,
                                'chr2': [chrmLabel] * num_entries,
                                'fragmentMid2': map(lambda i: mid_point_dict[index][i], np.nonzero(matrix)[1]),
                                'frag2': [1] * num_entries,
                                'score': matrix[np.nonzero(matrix)]
                               },
                               columns=['str1','chr1','fragmentMid1','frag1','str2','chr2','fragmentMid2','frag2','score'])
    Juicebox_df.to_csv('Juicebox/DKO_CD8_Juicebox_input.txt', sep='\t', index=None, header=None, mode='a')

for key in BD.transKeys:
    chr1, chr2 = key
    print chr1, chr2; sys.stdout.flush()
    index1 = BD._h5dict[str(key) + 'x']
    index2 = BD._h5dict[str(key) + 'y']
    values = BD._h5dict[str(key) + 'v']
    num_entries = len(index1)
    Juicebox_df = pd.DataFrame({'str1': [0] * num_entries,
                                'chr1': [BD.genome.idx2label[chr1]] * num_entries,
                                'fragmentMid1': map(lambda i: mid_point_dict[chr1][i], index1),
                                'frag1': [0] * num_entries,
                                'str2': [0] * num_entries,
                                'chr2': [BD.genome.idx2label[chr2]] * num_entries,
                                'fragmentMid2': map(lambda i: mid_point_dict[chr2][i], index2),
                                'frag2': [1] * num_entries,
                                'score': values
                                },
                                columns=['str1','chr1','fragmentMid1','frag1','str2','chr2','fragmentMid2','frag2','score'])
    Juicebox_df.to_csv('Juicebox/DKO_CD8_Juicebox_input.txt', sep='\t', index=None, header=None, mode='a')
