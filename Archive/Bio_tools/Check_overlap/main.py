#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  main.py
#  
#  Copyright 2017 Xiang <Xiang@LAPTOP-Q0TSHFKK>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Generral Public License as published by
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


import pandas as pd
import numpy as np

__RAW_DATA_PATH_DIR="/mnt/d/bioproject/Nanping_Check_overlap/atac-seq/raw_data/ignore_j_libraries/"
__Files_Name=["Spl2_CD4_young_j_peaks", "Thy2_CD4_young_j_peaks", "Spl2_CD8_young_j_peaks", "Thy2_CD8_young_j_peaks" ]
#__Files_Name=["test"]


def main():
	File_Type=".xlsx"

	for Fname in __Files_Name:
		File_Path=__RAW_DATA_PATH_DIR + Fname + File_Type
		print File_Path
		df1=pd.read_excel(File_Path, sheetname=0, usecols=[0,1,2,8], convert_float=True, header=0)	
		#df1=pd.read_csv(File_Path, usecols=[0,1,2], convert_float=True, header=0)
		#print (df1.iloc[:,:])
		df1=Manipulation_Excel(df1)
		Output(df1, Fname)
	
		print (df1.iloc[:,:])		

	return 0
	
def Manipulation_Excel(Input_File):
	#### Index of Drop rows
	drop_rows=[]
	Comparing_Line=3
	
	for rows in Input_File.index[:]:
		#print "Rows", rows
		if Input_File.iloc[rows, Comparing_Line] <= 2.99573227: #s (-LOG 0.05 = 2.99573227)
			drop_rows.append(rows)
	
	#Drop the column at position 1
	#df1_drop=df1.drop(df1.columns[[2]], axis=1)

	#### drop_row
	Input_File=Input_File.drop(Input_File.index[[drop_rows]])
	# Use `reset_index()` to reset the values
	Input_File=Input_File.reset_index(level=0, drop=True)
	#print Input_File
	
	return Input_File
	
def Output(Input_File, Output_Name):
	#Output
	File_Type=".bed"
	Out_Path=__RAW_DATA_PATH_DIR + Output_Name + "_Q_less_0.05" + File_Type
	# Create a Pandas Excel writer using XlsxWriter as the engine
	#Writer=pd.ExcelWriter(Out_Path, engine='xlsxwriter')
	# Convert the dataframe to an XlsxWriter Excel object.
	Input_File.to_csv(Out_Path, sep='	', mode='a', index=False, header=False)
	
	# Close the Pandas Excel writer and output the Excel file.
	#Input_File.save()


print "Start"
main()
print "End"




