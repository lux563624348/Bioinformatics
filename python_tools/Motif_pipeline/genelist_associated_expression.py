########################################################################
## 07/20/2018
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
########################################################################
# Usage 
#python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} \
#	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
########################################################################

import pandas as pd
import numpy as np
import sys
from optparse import OptionParser
import os
####################################################################################
## FUNCTIONS
### FUNCTION
def gene_associated_expression(genelist_path, file_name,expression_path, excelsheet_name):
    df = pd.read_excel(expression_path, sheet_name=excelsheet_name).drop_duplicates()
    df_genelist = pd.read_csv(genelist_path+'/'+file_name, sep='\t')
    df_genelist = df_genelist[['Gene Name','Distance to TSS']].rename(columns={'Gene Name':'gene_id'})
    df_genelist.merge(df, on='gene_id', how='inner').to_csv(genelist_path+'/Expression'+file_name, sep='\t', index=None)
    return 0

### For any input of /gene_exp.diff, return its up_DEGs genes with
### gene_id, cond1, cond2, log2(fold_change), p_value, plus a number of order.
### Parameters for DEGs:
FC_UP=2.0
q_value_less=0.05
FPKM_threshold=1.0
####################################################################################

def generate_Upregulated_Genes(Input_Path):
#### READ FILE
    df = pd.read_csv(Input_Path+'/gene_exp.diff', sep='\t', header=0, usecols={'test_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
#### Filter
    df=df[(df['status']=='OK') & (df['q_value']<=q_value_less) & (df['value_2']>=FPKM_threshold) & ( (df['log2(fold_change)'] >= np.log2(FC_UP)))]
#### Rearrange Columns
    df=df.rename(columns={'test_id':'gene_id','value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
#### Output
    return df.loc[:,['gene_id', df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)','p_value','q_value']]
####################################################################################
def generate_Downregulated_Genes(Input_Path):
#### READ FILE FROM CuffDiff Results
    df = pd.read_csv(Input_Path+'/gene_exp.diff', sep='\t', header=0, usecols={'test_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
#### Filter
    df=df[(df['status']=='OK') & (df['q_value']<=q_value_less) & (df['value_1']>=FPKM_threshold) & (df['log2(fold_change)']<= -np.log2(FC_UP))]   
#### Rearrange Columns
    df=df.rename(columns={'test_id':'gene_id','value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
#### Output
    return df.loc[:,['gene_id', df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)','p_value','q_value']]
####################################################################################
def diff_motif_MWU(motif_path, cuffdiff_path):
    import glob
    import fnmatch
    import scipy.stats as stats

    RAW_PATH=motif_path
    OUT_PATH=motif_path + '/Summary'
    DIR_CHECK_CREATE(OUT_PATH)
    FILE_TYPE='bed12.txt'

    df_up = generate_Upregulated_Genes (cuffdiff_path)
    df_down = generate_Downregulated_Genes(cuffdiff_path)

    df_output = pd.DataFrame(columns=['motif','pvalue'])
    i=0
    for path in glob.glob(RAW_PATH+'/*.'+FILE_TYPE):
        df_genelist = pd.read_csv(path,sep='\t')
        df_genelist = df_genelist[['Gene Name']].rename(columns={'Gene Name':'gene_id'}).drop_duplicates()
        motif_up_positive =  len(df_up.merge(df_genelist, on='gene_id',how='inner'))
        motif_up_negative = len(df_up) - motif_up_positive
        motif_down_positive = len(df_down.merge(df_genelist, on='gene_id',how='inner'))
        motif_down_negative = len(df_down) - motif_down_positive
        oddsratio, pvalue = stats.fisher_exact([[motif_up_positive,motif_down_positive],[motif_up_negative,motif_down_negative]])
        df_output.at[i,['motif','pvalue']] = path[-17:], pvalue
        i+=1
    df_output.to_csv(OUT_PATH+'/'+'Summary_on_Diff_Motif_associated_Expression.txt',sep='\t', index=None)
    return df_output
####################################################################################
def Add_common_header(df,common_header):
    name_list=[common_header]*len(df.columns)
    tuples = list(zip(name_list,df.columns))
    df.columns = pd.MultiIndex.from_tuples(tuples)
    return df
####################################################################################
def DIR_CHECK_CREATE(Input_Path):
	
    if (not os.path.isdir(Input_Path)):
        print ("New Dir Made is" + Input_Path)
        os.mkdir(Input_Path)
    else: 
        print ('Dir Exists')
    return 0
####################################################################################
### FUNCTION
### FUNCTIONS
def main(argv):
	desc="Associating genelist with of Epxression...."
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
		dest="cuffdiff", help="Path Contains Cuffdiff Results", metavar="<file>")
	parser.add_option("-m", "--motif", action="store", type="string",
		dest="motif_path", help="motif_path", metavar="<str>")
#	parser.add_option("-e", "--expression", action="store", type="string",
#		dest="expression_path", help="Full Path for Expression Excel File", metavar="<str>")
#	parser.add_option('-s', '--sheet_name', action="store", type="string",
#		dest="sheet_name", help="Sheet Name of Expression File", metavar="<str>")


	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	print " "
	print "Here is the Summary of your input."
	print "Input Path is %s" % opt.cuffdiff
	print "Input Motif Path is %s" % opt.motif_path
#	print "Input Expression Path is %s" % opt.expression_path
#	print "Input Expression Path is %s" % opt.sheet_name
	print " "
	
	diff_motif_MWU(opt.motif_path, opt.cuffdiff)
	#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
