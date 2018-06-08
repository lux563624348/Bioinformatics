## This file contains the basic function of read data and find the intersection and specific.
## Author: Xiang Li

###Sample
#### Input Package
import pandas as pd
import numpy as np
import os
import sys


INPUT_FOLDER_PATH=sys.argv[1]
INPUT_FILE_A=sys.argv[2]
INPUT_FILE_B=sys.argv[3]
FILE_TYPE=sys.argv[4]
Merge_Column=sys.argv[5]

def main():
	df_A=pd.read_csv(INPUT_FOLDER_PATH+INPUT_FILE_A+FILE_TYPE, sep='\t')
	df_B=pd.read_csv(INPUT_FOLDER_PATH+INPUT_FILE_B+FILE_TYPE, sep='\t')
	
	### Union A and B, use indicator to differentiate them.
	union_A_B=df_A.loc[:, [Merge_Column]].merge(df_B.loc[:,[Merge_Column]], how='outer', indicator=True)

	intersection_A_B=union_A_B[union_A_B['_merge']=='both']
	only_A=union_A_B[union_A_B['_merge']=='left_only']
	only_B=union_A_B[union_A_B['_merge']=='right_only']
	
	os.mkdir(INPUT_FOLDER_PATH+'/Results')
	#Output
	intersection_A_B.to_csv(INPUT_FOLDER_PATH+'/Results/Intersection_'+INPUT_FILE_A+'_'+INPUT_FILE_B+FILE_TYPE,index=None, columns=[Merge_Column])
	only_A.to_csv(INPUT_FOLDER_PATH+'/Results/Only_'+INPUT_FILE_A+FILE_TYPE,index=None, columns=[Merge_Column])
	only_B.to_csv(INPUT_FOLDER_PATH+'/Results/Only_'+INPUT_FILE_B+FILE_TYPE,index=None, columns=[Merge_Column])
	union_A_B.to_csv(INPUT_FOLDER_PATH+'/Results/Union_'+INPUT_FILE_A+'_'+INPUT_FILE_B+FILE_TYPE,index=None, columns=[Merge_Column])



print "Start"
main()
print "End"
