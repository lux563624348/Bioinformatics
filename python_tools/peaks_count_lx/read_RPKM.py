## Library: Vlad
## Author: Xiang Li

###Sample
#### Imput Package
import pandas as pd
import numpy as np
import sys
## If Input size = 0, using islandfiltered reads to normalize.



def main():
	INPUT_FILE_NAME = sys.argv[1]
	__EXE_PATH_DIR = sys.argv[2]
	__Input_Size = sys.argv[3]
	
	DATA_PATH = __EXE_PATH_DIR + '/' + INPUT_FILE_NAME + '.bed'
	
	print ("")
	print ("Calculate RPKM from bedtools intersect results.")
	print (DATA_PATH)
	print ("Input Library Size is " + __Input_Size)
	###format of bedtools gives,
	### 0   1   2    3          4
	###chr TSS TES gene_id   Read_count
	df_tem = pd.read_csv(DATA_PATH, sep='\t', header=-1)
	df_tem = df_tem.rename(columns={0: "chr", 1: "TSS", 2: "TES", 3: "gene_id", 4: "#_reads"})
	df_tem['len'] = np.abs(df_tem.loc[:,'TSS'] - df_tem.loc[:,'TES'])
	__Input_Size=float(__Input_Size)
	if(__Input_Size==0):
		print "Input_Size = 0, using islandreads to normalize"
		number_of_total_reads = float(df_tem[['#_reads']].sum(axis=0))
	else:
		number_of_total_reads = __Input_Size
	RPKM=df_tem.loc[:,'#_reads']*(10.0**9)/((df_tem.loc[:,'len'])*number_of_total_reads)

	df = pd.DataFrame(columns=['gene_id','RPKM'])
	df['gene_id'] = df_tem['gene_id']
	df.loc[:,'RPKM'] = RPKM.values
	df = df.set_index('gene_id')
	
	OUT_PATH= __EXE_PATH_DIR + '/RPKM_' + INPUT_FILE_NAME + '.csv'
	
	df.to_csv(OUT_PATH, sep='\t')
	

print "Start"
main()
print "End"
