import os
import numpy as np
from mirnylib import h5dict, genome
from hiclib import mapping, fragmentHiC, binnedData

genome_db = genome.Genome('/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jun2016/exp/iterative_mapping/genome.fa')

mapped_reads = h5dict.h5dict('DKO_CD8_p1_mapped_reads.hdf5')

mapping.parse_sam(
    sam_basename1='DKO_CD8/DKO_CD8_R1_p1.bam',
    sam_basename2='DKO_CD8/DKO_CD8_R2_p1.bam',
    out_dict=mapped_reads,
    genome_db=genome_db)
