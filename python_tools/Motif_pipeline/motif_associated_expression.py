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
import sys
import numpy as np
from optparse import OptionParser
import os
####################################################################################
## FUNCTIONS
FC_UP=2.0
q_value_less=0.05
FPKM_threshold=1.0
####################################################################################
### FUNCTION
def associate_motif_expression_txt(motif_gene_association_input, expression_input):
    df_motif_gene = pd.read_csv(motif_gene_association_input, sep='\t')
    df_gene_exp = pd.read_csv(expression_input, sep='\t')
    df_gene_exp = df_gene_exp.rename(columns={'test_id':'Gene Name'}).drop(['gene_id','gene'],axis=1)
    df_merge = df_motif_gene.merge(df_gene_exp,on='Gene Name', how='inner')
    gene_length = df_merge['locus'].str.split(pat=':',expand=True)[1].str.split(pat='-',expand=True).astype(int).apply(lambda x: x[1] - x[0], axis=1)
    df_merge.insert(5,'Distance to TES', df_merge['Distance to TSS'] - gene_length)
    df_merge.to_csv(motif_gene_association_input[:-4]+'_expression.txt',sep='\t', index=False)
    print "Motif Associated with Exp Finished!"
    return 0

def generate_Upregulated_Genes(df):
#### Filter
    df=df[(df['status']=='OK') & (df['q_value']<=q_value_less) & (df['value_2']>=FPKM_threshold) & ( (df['log2(fold_change)'] >= np.log2(FC_UP)))]
#### Output
    return df
####################################################################################
def generate_Downregulated_Genes(df):
#### Filter
    df=df[(df['status']=='OK') & (df['q_value']<=q_value_less) & (df['value_1']>=FPKM_threshold) & (df['log2(fold_change)']<= -np.log2(FC_UP))]   
#### Output
    return df
   
def associate_motif_expression_excel(motif_gene_association_folder, expression_input):
    import glob
    import datetime
    from scipy import stats

    os.chdir(motif_gene_association_folder)

    INPUT_LIST = glob.glob('*.bed.gene_association.txt')
    writer = pd.ExcelWriter(motif_gene_association_folder+'/Motif_Associated_Genes_Summary_'+str(datetime.date.today())+ '.xlsx', engine='xlsxwriter')
    df_gene_exp = pd.read_csv(expression_input, sep='\t').rename(columns={'test_id':'Gene Name'}).drop(['gene_id','gene','test_stat'],axis=1)

    
    summary_df = pd.DataFrame(columns=['Motif_Name','wilcoxon_rank_test_pvalue','#_Up','#_Down','#_All','(Up-Down)/All']) ### wilcoxon signed rank test 
    i=0
    for name in INPUT_LIST[:]:
        df_motif_gene = pd.read_csv(name, sep='\t')
        df_merge = df_motif_gene.merge(df_gene_exp,on='Gene Name', how='inner')
        ## Gene Length Calculation
        gene_length = df_merge['locus'].str.split(pat=':',expand=True)[1].str.split(pat='-',expand=True).astype(int).apply(lambda x: x[1] - x[0], axis=1)
        df_merge.insert(5,'Distance to TES', df_merge['Distance to TSS'] - gene_length)
        ## Save a summary of motif associated with Expression
        Cond1_Name = df_merge['sample_1'][0]
        Cond2_Name = df_merge['sample_2'][0]
        df_merge.rename(columns={'value_1':Cond1_Name, 'value_2':Cond2_Name}).drop(['sample_1','sample_2'],
                axis=1).to_excel(writer, sheet_name=name[0:7], index=None)

        ### https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ks_2samp.html
        ### KS 2 sample test
        ### https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html
        #print "pvalue:",stats.ks_2samp(df_merge['value_1'], df_merge['value_2'])[1] ## wilcoxon signed rank test 
        num_up = len(generate_Upregulated_Genes(df_merge)['Gene Name'].unique())
        num_down = len(generate_Downregulated_Genes(df_merge)['Gene Name'].unique())
        num_all = len(df_merge['Gene Name'].unique())
        up_ratio = 1.0*(num_up -num_down  )/num_all
        summary_df.loc[i,['Motif_Name','wilcoxon_rank_test_pvalue','#_Up','#_Down','#_All','(Up-Down)/All']] =\
                       [name[0:7], stats.wilcoxon(df_merge['value_1'], df_merge['value_2'])[1],num_up, num_down, num_all, up_ratio]
        i+=1
    #num_up = len(generate_Upregulated_Genes(df_gene_exp)['Gene Name'])
    #num_down = len(generate_Downregulated_Genes(df_gene_exp)['Gene Name'])
    #num_all = len(df_gene_exp[df_gene_exp['status']=='OK']['Gene Name'].unique())  ##Choose Gene Expression satisfied statistical significance
    #up_ratio =  1.0*(num_up-num_down)  / num_all
    #summary_df.loc[i,['Motif_Name','wilcoxon_rank_test_pvalue','#_Up','#_Down','#_All','(Up-Down)/All']] =\
    #['Overall', stats.wilcoxon(df_gene_exp['value_1'], df_gene_exp['value_2'])[1],num_up, num_down, num_all, up_ratio]
    summary_df.sort_values(by=['(Up-Down)/All'],ascending=False).to_excel(writer, sheet_name='Motif Ranking Test', index=None)
    writer.save()
    
    return 0
####################################################################################
### FUNCTION
### FUNCTIONS
def main(argv):
	desc="Associating Motif & Gene with of Epxression...."
	parser = OptionParser(description=desc)
	parser.add_option("-e", "--expression", action="store", type="string",
		dest="gene_exp", help="Path of Gene Expression from Cuffdiff", metavar="<file>")
	parser.add_option("-m", "--motif", action="store", type="string",
		dest="motif_gene", help="motif with genes", metavar="<dir>")


	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	print " "
	print "Here is the Summary of your input."
	print "Input Gene Expression Path is %s" % opt.gene_exp
	print "Input Motif & Gene Path is %s" % opt.motif_gene
	print " "
	
	associate_motif_expression_excel(opt.motif_gene, opt.gene_exp)




if __name__ == "__main__":
	main(sys.argv)
