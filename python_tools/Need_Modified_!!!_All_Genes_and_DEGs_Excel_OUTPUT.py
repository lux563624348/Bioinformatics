#!/usr/bin/env python
### From CuffDiff Results, generate DEGs and All_Genes in one excel file.
#  

## This file contains the basic function of read data and clean filter.
## Author: Xiang Li

###Sample
#### Imput Package
import pandas as pd
import numpy as np
import os



def main(args):
	#### Generate INPUT PATH
	os.makedirs(os.getcwd()+'/genelist/')
	PATH_FOLDER=os.getcwd()+ '/CuffDiff_Results/'
	OUT_FOLDER=os.getcwd()+'/genelist/'


	INPUT_LIST=os.listdir(PATH_FOLDER)
	writer = pd.ExcelWriter(OUT_FOLDER+'All_genes_DEGs.xlsx', engine='xlsxwriter')

	i=0
	for input_name in INPUT_LIST:
		INPUT_PATH = PATH_FOLDER+input_name+'/gene_exp.diff'
		if (i==0):
			df_all=generate_All_Genes(INPUT_PATH,i)
			print(df_all.shape)
			i+=1
			continue
		df_all = df_all.merge(generate_All_Genes(INPUT_PATH,i), on='gene_id', how='inner') 
		print ('# of Up:' )
		print(df_all.shape)
		i+=1
		
	df_all=df_all.set_index('gene_id')
	df_all.to_excel(writer, sheet_name='All_Genes')
		
	PATH_FOLDER=os.getcwd()+ '/CuffDiff_Results/'
	OUT_FOLDER=os.getcwd()+'/genelist/'

	for input_name in INPUT_LIST:
		INPUT_PATH = PATH_FOLDER+input_name+'/gene_exp.diff'
		df_up = generate_Upregulated_Genes(INPUT_PATH)
		df_up.to_excel( writer, sheet_name=input_name)

		df_down = generate_Downregulated_Genes(INPUT_PATH)
		df_down.to_excel( writer, sheet_name=input_name)
	writer.save()
		
		
		
		
		
		
		
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main(sys.argv))



def generate_All_Genes(Input_Path, number):
    #### READ FILE
    df = pd.read_csv(Input_Path, sep='\t', header=0, usecols={'gene_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
    #### Filter
# None
    #### Rearrange Columns
    
    df=df.rename(columns={'value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
    #df=df.set_index('gene_id')
    #### Output
    
    return df.loc[:,['gene_id',df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)','p_value','q_value', number]].fillna('')
    
def generate_Upregulated_Genes(Input_Path):
    #### READ FILE
    df = pd.read_csv(Input_Path, sep='\t', header=0, usecols={'gene_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
    #### Filter
    
    df=df[
    (df['status']=='OK') & (df['q_value']<0.05) & (df['value_2']>=1.0)\
    & (df['log2(fold_change)'] > np.log2(1.5))
    ]
    
    #### Rearrange Columns
    
    df=df.rename(columns={'value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
    df=df.set_index('gene_id')
    #### Output
    
    return df.loc[:,[df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)','p_value','q_value']]

def generate_Downregulated_Genes(Input_Path):
    #### READ FILE FROM CuffDiff Results
    df = pd.read_csv(Input_Path, sep='\t', header=0, usecols={'gene_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
    #### Filter
    df=df[
    (df['status']=='OK') & (df['q_value']<0.05) & (df['value_1']>=1.0)\
    & (df['log2(fold_change)']< np.log2( 2.0/3.0))
    ]
    
    #### Rearrange Columns
    
    df=df.rename(columns={'value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
    df=df.set_index('gene_id')
    #### Output
    return df.loc[:,[df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)','p_value','q_value']]
