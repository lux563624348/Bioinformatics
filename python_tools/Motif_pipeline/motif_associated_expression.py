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
from optparse import OptionParser
import os
####################################################################################
## FUNCTIONS
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
    
def associate_motif_expression_excel(motif_gene_association_folder, expression_input):
	import glob
	import datetime
	os.chdir(motif_gene_association_folder)
	INPUT_LIST = glob.glob('*.bed.gene_association.txt')
	
	writer = pd.ExcelWriter(motif_gene_association_folder+'/Motif_Associated_Genes_Summary_'+str(datetime.date.today())+ '.xlsx', engine='xlsxwriter')
	df_gene_exp = pd.read_csv(expression_input, sep='\t').rename(columns={'test_id':'Gene Name'}).drop(['gene_id','gene'],axis=1)
	
	for name in INPUT_LIST[:]:
		df_motif_gene = pd.read_csv(name, sep='\t')
		df_merge = df_motif_gene.merge(df_gene_exp,on='Gene Name', how='inner')
		gene_length = df_merge['locus'].str.split(pat=':',expand=True)[1].str.split(pat='-',expand=True).astype(int).apply(lambda x: x[1] - x[0], axis=1)
		df_merge.insert(5,'Distance to TES', df_merge['Distance to TSS'] - gene_length)
		df_merge.to_excel(writer, sheet_name=name[0:7], index=None)
	writer.save()
	return "Motif Associated with Exp"
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
