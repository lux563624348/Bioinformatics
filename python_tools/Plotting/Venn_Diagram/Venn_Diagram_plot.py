#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  
#  Copyright 2017 Xiang <Xiang@LAPTOP-Q0TSHFKK>
#  
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import sys

INPUT_FILE=sys.argv[1]
ConA=sys.argv[2]
ConB=sys.argv[3]
ConC=sys.argv[4]

def main():
	df1=pd.read_csv(INPUT_FILE, delimiter='	', header=-1, \
	names= ['Chr', 'TSS','TES', ConA, ConB, ConC])
	
	A_Sum=df1[ConA].sum()
	B_Sum=df1[ConB].sum()
	C_Sum=df1[ConC].sum()
	
	df1.loc[df1.index[df1[ConA]>1], [ConA]]=1
	df1.loc[df1.index[df1[ConB]>1], [ConB]]=1
	df1.loc[df1.index[df1[ConC]>1], [ConC]]=1
	
	df1['ABC']=df1[ConA]+df1[ConB]+df1[ConC]
	df1['AB']=df1[ConA]+df1[ConB]
	df1['BC']=df1[ConB]+df1[ConC]
	df1['AC']=df1[ConA]+df1[ConC]
	
	center=df1['ABC'].value_counts()[3]
	AB_com=df1['AB'].value_counts()[2]-center
	AC_com=df1['AC'].value_counts()[2]-center
	BC_com=df1['BC'].value_counts()[2]-center


	A_solo=A_Sum-AB_com-AC_com-center
	B_solo=B_Sum-AB_com-BC_com-center
	C_solo=C_Sum-AC_com-BC_com-center


	v=venn3(subsets = (A_solo, B_solo, AB_com, C_solo, AC_com, BC_com,center), \
	set_labels = (ConA, ConB, ConC))
	#c=venn3_circles(subsets = (10, 8, 22, 6,9,4,2), linestyle='dashed', linewidth=1, color="grey")
	plt.title("Venn Diagram of 3 sets")
	plt.savefig('Venn.png')

	
	#print (Con1)
	#print (Con2)
	#print (Con3)
	
print "Start"
main()
print "End"
