#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  read_count_from_bedpeaks.py
#  
#  Copyright 2017 Xiang <Xiang@LAPTOP-Q0TSHFKK>
#  

import sys
import numpy as np
import pandas as pd
#from optparse import OptionParser

__RAW_DATA_PATH_DIR = sys.argv[1]
__EXE_PATH = sys.argv[2]
Output_name = sys.argv[3]

def main():
	ref_gene_data = pd.read_csv(sys.argv[1]+ '/OUTPUT/genes.bed', delimiter='	',header=-1, names= ['Chr', 'Up_4k' ,'TES', 'Gene_name','N_Peaks']) #, nrows=100
	ref_peak_data = pd.read_csv(sys.argv[1]+ '/OUTPUT/peaks_counting.bed', delimiter='	', header=-1, names= ['Chr', 'TSS' ,'TES','N_reads', 'p_value', 'n_overlap']) #, nrows=5321)

	Total_reads = ref_peak_data.N_reads.sum()
	Total_Peaks = ref_gene_data.N_Peaks.sum()
	
	
	ref_gene_data = ref_gene_data[ ref_gene_data.N_Peaks != 0 ]
	ref_gene_data = ref_gene_data.reset_index(level=None)
	
	ref_peak_data = ref_peak_data[ ref_peak_data.n_overlap != 0 ]
	
	#Total_Peaks = ref_gene_data.N_Peaks.sum()
	#Total_overlap = ref_peak_data.n_overlap.sum()
	
	Start_gene_order=0
	Reads_count = pd.Series(np.zeros(len(ref_gene_data)))
	Peaks_count = pd.Series(np.zeros(len(ref_gene_data)))
	RPKM_count = pd.Series(np.zeros(len(ref_gene_data)))
	gene_length= pd.Series(np.zeros(len(ref_gene_data)))
	
	#print (len(ref_peak_data))
	
	
	for i in range(len(ref_peak_data)):
		while ( Peaks_count.iloc[Start_gene_order] == ref_gene_data.N_Peaks.values [Start_gene_order] ):
			Start_gene_order += 1
		buff = 0
		for j in range(ref_peak_data.n_overlap.values[i]):
			if ( Peaks_count.iloc[Start_gene_order + j + buff ] < ref_gene_data.N_Peaks.values [Start_gene_order + j + buff] ):
				Reads_count.iloc[ Start_gene_order + j + buff ] += ref_peak_data.N_reads.values[i]
				Peaks_count.iloc[ Start_gene_order + j + buff ] += 1
			else:
				while ( Peaks_count.iloc[Start_gene_order + j + buff ] == ref_gene_data.N_Peaks.values [Start_gene_order + j + buff] ):
					buff +=1
				else:
					Reads_count.iloc[ Start_gene_order + j + buff ] += ref_peak_data.N_reads.values[i]
					Peaks_count.iloc[ Start_gene_order + j + buff ] += 1
					
					
	for i in range(len(ref_gene_data)):
		gene_length.iloc[i] = abs( ref_gene_data.Up_4k.values[i] - ref_gene_data.TES.values[i] )
		RPKM_count.iloc[i] = Reads_count.iloc[i]*10.0**9 /(gene_length.iloc[i] * Total_reads)

	
	Out_put = pd.DataFrame( { '1.NAME': ref_gene_data.Gene_name, '2.Length': gene_length, '3.Num_Peaks': ref_gene_data.N_Peaks, '4.Num_Reads': Reads_count, '5.RPKM': RPKM_count }) 
	Output(Out_put, Output_name)
	#print gene_length
	
	#Check Point
	#print ('4930594C11Rik: 108, Peaks: 3', Reads_count.iloc[76], Peaks_count.iloc[76])
	#print ('Opn3: 196, Peaks: 6', Reads_count.iloc[533], Peaks_count.iloc[533])
	#print ('Chml: 89, Peaks: 2', Reads_count.iloc[534], Peaks_count.iloc[534])
	#print ('Wdr64: 998, Peaks: 20', Reads_count.iloc[535], Peaks_count.iloc[535])
	#print ('Raet1a: 3276, Peaks: 42', Reads_count.iloc[664], Peaks_count.iloc[664])
	#print ('Raet1e: 38, Peaks: 2', Reads_count.iloc[665], Peaks_count.iloc[665])
	#print ('Raet1b: 3106, Peaks: 39', Reads_count.iloc[666], Peaks_count.iloc[666])
	#print ('Raet1c: 3106, Peaks: 39', Reads_count.iloc[667], Peaks_count.iloc[667])
	#print ('H60b: 1079, Peaks: 13', Reads_count.iloc[668], Peaks_count.iloc[668])
	#print ('i',i, 'ref_peak_data.N_reads' ,ref_peak_data.N_reads.values[i])
	
	return 0

	
def Output(Input_File, Output_Name):
	File_Type=".bed"
	Out_Path= __EXE_PATH + "/RPKM_" + Output_Name + File_Type
	Input_File.to_csv(Out_Path, sep='	', mode='w', index=False, header=True)
	return 0

print "Start"
main()
print "End"








