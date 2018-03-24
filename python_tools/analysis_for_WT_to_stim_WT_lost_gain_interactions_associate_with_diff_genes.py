#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  analysis_for_WT_to_stim_WT_lost_gain_interactions_associate_with_diff_genes.py
#  
#  Copyright 2018 Xiang <Xiang@LAPTOP-Q0TSHFKK>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


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
import HTSeq




def main(args):
	genome_db = genome.Genome('/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM9/genome.fa')
	
	
	up_genes_dict = {}
	up_genes_bed_file = HTSeq.BED_Reader('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis/processed/HiC_around_Tcf1_peaks/diff_genes/HP_WT_up_gene_promoter_iv.bed')
	for alt in up_genes_bed_file:
		up_genes_dict[alt.name] = alt.iv
	
	down_genes_dict = {}
	down_genes_bed_file = HTSeq.BED_Reader('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis/processed/HiC_around_Tcf1_peaks/diff_genes/HP_WT_down_gene_promoter_iv.bed')
	for alt in down_genes_bed_file:
		down_genes_dict[alt.name] = alt.iv
	
	
	
	
	global idx
	idx=pd.IndexSlice
	
	Con1_with_cutoff()
	Con2_with_cutoff()
	Con1_all()
	Con2_all()
	
	stim_WT_sig_df.shape
	stim_WT_all_df.shape
	stim_DKO_sig_df.head()
	
	
	Tcf1_CD8_interactions_df = Tcf1_WT_interactions_df.merge(Tcf1_stim_DKO_interactions_df, on=['Tcf1_binding_iv', 'chrom', 'Tcf1_binding_bin_start', 'interaction_site_bin_start'], how='left')
	Tcf1_CD8_interactions_df['stim_DKO_CD8_score'] = Tcf1_CD8_interactions_df.apply(lambda row: retrieve_interaction_score_from_stim_DKO(row, 'score'), axis=1)
	Tcf1_CD8_interactions_df['stim_DKO_CD8_balanced_count'] = Tcf1_CD8_interactions_df.apply(lambda row: retrieve_interaction_score_from_stim_DKO(row, 'balanced_count'), axis=1)
	
	return 0

def iv_to_str(iv):
		return iv.chrom + ':' + str(iv.start) + '-' + str(iv.end)  

def retrieve_interaction_score_from_stim_DKO(row, colname):
		start1 = min(row['Tcf1_binding_bin_start'], row['interaction_site_bin_start'])
		start2 = max(row['Tcf1_binding_bin_start'], row['interaction_site_bin_start'])
		try:
			return stim_DKO_all_df.loc[idx[row['chrom'], start1, start2], colname]
		except:
			return 0


def Con1_with_cutoff():
	os.chdir('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis')

	cutoff = 10e-5
	sig_df = pd.read_csv('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/fit-hi-c/stim_WT_CD8.spline_pass2.significances.txt.gz', header=0, sep='\t')
	sig_df = sig_df[sig_df['q-value'] < cutoff]

	mycooler = cooler.Cooler('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/stim_WT_CD8-res-10K.cool')
	bin_df = mycooler.bins()[:]
	bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)

	bin_start_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_start_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'start']].values)

	bin_end_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_end_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'end']].values)

	sig_df['start1'] = sig_df.apply(lambda row: bin_start_dict[row['chr1']][row['fragmentMid1']], axis=1)
	sig_df['start2'] = sig_df.apply(lambda row: bin_start_dict[row['chr2']][row['fragmentMid2']], axis=1)

	sig_df['name'] = sig_df.index
	sig_df['strand'] = '.'
	sig_df['score'] = -np.log10(sig_df['p-value'])
	sig_df.loc[sig_df.score == np.inf, 'score'] = 1000

	chr_matrix_list = []
	for chrom in mycooler.chromnames:
		chr_matrix = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chrom).fillna(0)
		chr_matrix = chr_matrix[(abs(chr_matrix.start1 - chr_matrix.start2) >= 20000) & (abs(chr_matrix.start1 - chr_matrix.start2) <= 2000000)]
		chr_matrix_list.append(chr_matrix)

	matrix_df = pd.concat(chr_matrix_list).set_index(['chrom1', 'start1', 'start2'])
	sig_df['balanced_count'] = sig_df.apply(lambda row: matrix_df.loc[idx[row['chr1'], row['start1'], row['start2']], 'balanced'], axis=1)

	global stim_WT_sig_df
	stim_WT_sig_df = sig_df

	stim_WT_sig_df.set_index(['chr1', 'start1', 'start2'], inplace=True)
	
def Con2_with_cutoff():
	os.chdir('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis')

	cutoff = 10e-5
	sig_df = pd.read_csv('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/fit-hi-c/stim_DKO_CD8.spline_pass2.significances.txt.gz', header=0, sep='\t')
	sig_df = sig_df[sig_df['q-value'] < cutoff]

	mycooler = cooler.Cooler('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/stim_DKO_CD8-res-10K.cool')
	bin_df = mycooler.bins()[:]
	bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)

	bin_start_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_start_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'start']].values)

	bin_end_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_end_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'end']].values)

	sig_df['start1'] = sig_df.apply(lambda row: bin_start_dict[row['chr1']][row['fragmentMid1']], axis=1)
	sig_df['start2'] = sig_df.apply(lambda row: bin_start_dict[row['chr2']][row['fragmentMid2']], axis=1)

	sig_df['name'] = sig_df.index
	sig_df['strand'] = '.'
	sig_df['score'] = -np.log10(sig_df['p-value'])
	sig_df.loc[sig_df.score == np.inf, 'score'] = 1000

	chr_matrix_list = []
	for chrom in mycooler.chromnames:
		chr_matrix = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chrom).fillna(0)
		chr_matrix = chr_matrix[(abs(chr_matrix.start1 - chr_matrix.start2) >= 20000) & (abs(chr_matrix.start1 - chr_matrix.start2) <= 2000000)]
		chr_matrix_list.append(chr_matrix)

	matrix_df = pd.concat(chr_matrix_list).set_index(['chrom1', 'start1', 'start2'])
	sig_df['balanced_count'] = sig_df.apply(lambda row: matrix_df.loc[idx[row['chr1'], row['start1'], row['start2']], 'balanced'], axis=1)

	global stim_DKO_sig_df 
	stim_DKO_sig_df = sig_df

	stim_DKO_sig_df.set_index(['chr1', 'start1', 'start2'], inplace=True)

def Con1_all():
	os.chdir('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis')

	sig_df = pd.read_csv('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/fit-hi-c/stim_WT_CD8.spline_pass2.significances.txt.gz', header=0, sep='\t')

	mycooler = cooler.Cooler('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/stim_WT_CD8-res-10K.cool')
	bin_df = mycooler.bins()[:]
	bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)

	bin_start_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_start_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'start']].values)

	bin_end_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_end_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'end']].values)

	sig_df['start1'] = sig_df.apply(lambda row: bin_start_dict[row['chr1']][row['fragmentMid1']], axis=1)
	sig_df['start2'] = sig_df.apply(lambda row: bin_start_dict[row['chr2']][row['fragmentMid2']], axis=1)

	sig_df['name'] = sig_df.index
	sig_df['strand'] = '.'
	sig_df['score'] = -np.log10(sig_df['p-value'])
	sig_df.loc[sig_df.score == np.inf, 'score'] = 1000

	chr_matrix_list = []
	for chrom in mycooler.chromnames:
		chr_matrix = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chrom).fillna(0)
		chr_matrix = chr_matrix[(abs(chr_matrix.start1 - chr_matrix.start2) >= 20000) & (abs(chr_matrix.start1 - chr_matrix.start2) <= 2000000)]
		chr_matrix_list.append(chr_matrix)

	matrix_df = pd.concat(chr_matrix_list).set_index(['chrom1', 'start1', 'start2'])
	sig_df['balanced_count'] = sig_df.apply(lambda row: matrix_df.loc[idx[row['chr1'], row['start1'], row['start2']], 'balanced'], axis=1)

	global stim_WT_all_df
	stim_WT_all_df = sig_df

	stim_WT_all_df.set_index(['chr1', 'start1', 'start2'], inplace=True)
	
def Con2_all():
	os.chdir('/home/lxiang/cloud_research/PengGroup/XLi/Data/Haihui/Tcf1/HiC-seq/Analysis')

	sig_df = pd.read_csv('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/fit-hi-c/stim_DKO_CD8.spline_pass2.significances.txt.gz', header=0, sep='\t')

	mycooler = cooler.Cooler('/home/lxiang/cloud_research/PengGroup/ZZeng/Data/Haihui/Tcf1/HiC-seq/Jul2017/iterative_mapping/stim_DKO_CD8-res-10K.cool')
	bin_df = mycooler.bins()[:]
	bin_df['mid'] = ((bin_df['start'] + bin_df['end']) / 2).astype(int)

	bin_start_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_start_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'start']].values)

	bin_end_dict = collections.defaultdict(lambda: {})
	for chrom in mycooler.chromnames:
		bin_end_dict[chrom] = dict(bin_df.loc[bin_df.chrom == chrom, ['mid', 'end']].values)

	sig_df['start1'] = sig_df.apply(lambda row: bin_start_dict[row['chr1']][row['fragmentMid1']], axis=1)
	sig_df['start2'] = sig_df.apply(lambda row: bin_start_dict[row['chr2']][row['fragmentMid2']], axis=1)

	sig_df['name'] = sig_df.index
	sig_df['strand'] = '.'
	sig_df['score'] = -np.log10(sig_df['p-value'])
	sig_df.loc[sig_df.score == np.inf, 'score'] = 1000

	chr_matrix_list = []
	for chrom in mycooler.chromnames:
		chr_matrix = mycooler.matrix(balance=True, as_pixels=True, join=True).fetch(chrom).fillna(0)
		chr_matrix = chr_matrix[(abs(chr_matrix.start1 - chr_matrix.start2) >= 20000) & (abs(chr_matrix.start1 - chr_matrix.start2) <= 2000000)]
		chr_matrix_list.append(chr_matrix)

	matrix_df = pd.concat(chr_matrix_list).set_index(['chrom1', 'start1', 'start2'])
	sig_df['balanced_count'] = sig_df.apply(lambda row: matrix_df.loc[idx[row['chr1'], row['start1'], row['start2']], 'balanced'], axis=1)

	global stim_DKO_all_df
	stim_DKO_all_df = sig_df

	stim_DKO_all_df.set_index(['chr1', 'start1', 'start2'], inplace=True)
	
def Con1_interactions():
	Tcf1_stim_WT_interactions_df_list = []

	score_cutoff = 0
	cutoff_sig_df = stim_WT_sig_df[stim_WT_sig_df.score >= score_cutoff]
	cutoff_sig_1_df = cutoff_sig_df.set_index(['chr1', 'start1', 'start2'])
	cutoff_sig_2_df = cutoff_sig_df.set_index(['chr1', 'start2', 'start1'])

	for peak, peak_iv in peaks_dict.iteritems():
		peak_mid = peak_iv.start
		peak_bin_start = peak_mid / 10000 * 10000
		chrom = peak_iv.chrom
		
		try:
			peak_cutoff_sig_2_df = cutoff_sig_2_df.loc[idx[chrom, peak_bin_start],['score', 'balanced_count']]
			peak_cutoff_sig_2_df['interaction_site_bin_start'] = peak_cutoff_sig_2_df.index
			peak_cutoff_sig_2_df['Tcf1_binding_bin_start'] = peak_bin_start
			peak_cutoff_sig_2_df['chrom'] = chrom
			peak_cutoff_sig_2_df['Tcf1_binding_iv'] = peak
		except:
			peak_cutoff_sig_2_df = None
		
		try:
			peak_cutoff_sig_1_df = cutoff_sig_1_df.loc[idx[chrom, peak_bin_start],['score', 'balanced_count']]
			peak_cutoff_sig_1_df['interaction_site_bin_start'] = peak_cutoff_sig_1_df.index
			peak_cutoff_sig_1_df['Tcf1_binding_bin_start'] = peak_bin_start
			peak_cutoff_sig_1_df['chrom'] = chrom
			peak_cutoff_sig_1_df['Tcf1_binding_iv'] = peak
		except:
			peak_cutoff_sig_1_df = None
		
		if peak_cutoff_sig_1_df is not None or peak_cutoff_sig_2_df is not None:
			Tcf1_stim_WT_interactions_df_list.append(pd.concat([peak_cutoff_sig_1_df, peak_cutoff_sig_2_df])) 
	Tcf1_stim_WT_interactions_df = pd.concat(Tcf1_stim_WT_interactions_df_list).loc[:,['Tcf1_binding_iv', 'chrom', 'Tcf1_binding_bin_start', 'interaction_site_bin_start', 'score', 'balanced_count']]
	Tcf1_stim_WT_interactions_df.rename(columns={'score': 'stim_WT_CD8_score', 'balanced_count': 'WT_CD8_balanced_count'}, inplace=True)
	
def 
	
if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

print "Start"
main()
print "End"
