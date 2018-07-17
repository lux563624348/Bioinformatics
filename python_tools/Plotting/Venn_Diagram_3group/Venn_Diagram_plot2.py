#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  
#  Copyright 2017 Xiang <Xiang@LAPTOP-Q0TSHFKK>
#
'''
For example:
Chr	 start	 end
File A: Chr1 	1	20
File B: Chr1	15	30
File C: Chr1	45	60
> Union> File: union.bed
Chr	 start	 end
Chr1 	1	30
Chr1	45	60

Overlap union.bed with A, B and C one by one.
				‘A’	‘B’	‘C’
Chr1 	1	30	1	1	0  
Chr1	45	60	0	0	0

Because chr1	1	30 has overlap both in ‘A’ and ‘B’, then we says that this peak is common in A and B.
'''

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys

INPUT_FILE=sys.argv[1]
ConA=sys.argv[2]
ConB=sys.argv[3]

def main():
	df1=pd.read_csv(INPUT_FILE, delimiter='	', header=-1, \
	names= ['Chr', 'TSS','TES', ConA, ConB])
	
	A_Sum=df1[ConA].sum()
	B_Sum=df1[ConB].sum()
	
	df1.loc[df1.index[df1[ConA]>1], [ConA]]=1
	df1.loc[df1.index[df1[ConB]>1], [ConB]]=1
	
	df1['AB']=df1[ConA]+df1[ConB]

	AB_com=df1['AB'].value_counts()[2]

	A_solo=A_Sum-AB_com
	B_solo=B_Sum-AB_com



	v=venn2(subsets = (A_solo, B_solo, AB_com), \
	set_labels = (ConA, ConB))
	#c=venn3_circles(subsets = (10, 8, 22, 6,9,4,2), linestyle='dashed', linewidth=1, color="grey")
	plt.title("Venn Diagram of 2 sets")
	plt.savefig('Venn2.png')

	
print "Start"
main()
print "End"
