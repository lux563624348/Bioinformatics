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
import os
from optparse import OptionParser

####################################################################################
## FUNCTIONS
### FUNCTION
def gene_associated_expression(genelist_path, file_name,expression_path, excelsheet_name):
    df = pd.read_excel(expression_path, sheet_name=excelsheet_name).drop_duplicates()
    df_genelist = pd.read_csv(genelist_path+'/'+file_name, sep='\t')
    df_genelist = df_genelist[['Gene Name','Distance to TSS']].rename(columns={'Gene Name':'gene_id'})
    df_genelist.merge(df, on='gene_id', how='inner').to_csv('Expression'+file_name, sep='\t', index=None)
    return 0

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
	desc="After CuffDiff, summary of Epxression and PCA heatmap and pearson correlation analysis."
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="cuffdiff_path", help="Path Contains Folder CuffDiff_Results", metavar="<file>")
	parser.add_option("-n", "--project_name", action="store", type="string",
		dest="project_name", help="Project Name for Plot", metavar="<str>")
	parser.add_option("-l", "--list_replicates_orders", action="store", type="int",
		dest="replicates_orders", help="Number of Replicates per lib. Such as 132 means one rep for first lib, 3 replicates for second and 2 for third.", metavar="<str>")
	parser.add_option('-f', '--fold_change', action="store", type="float",dest="fold_change", 
		help="Fold Change to Define Differential Expression Genes", metavar="<float>")
	parser.add_option('-q', '--q_value', action="store", type="float",dest="q_value", 
		help="Required q_value", metavar="<float>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 5:
		parser.print_help()
		sys.exit(1)
	
	global FC_UP,q_value_less
	FC_UP = opt.fold_change
	q_value_less = opt.q_value
	replicates_orders_list =  [int(x) for x in str(opt.replicates_orders)]
	
	print " "
	print "Here is the Summary of your input."
	print "Input Path Contains CuffDiff_Results : %s" % opt.cuffdiff_path
	print "Project Name for Plot: %s" % opt.project_name
	print "Number of replicates per lib: %s" % replicates_orders_list
	print "old Change to Define Differential Expression Genes: %f" % opt.fold_change
	print "Required q_value: %f" % opt.q_value
	print "End of Summary."
	print " "
	
	### Two steps
	df_DEGs, df_ALL_Genes = Summary_Expression( opt.cuffdiff_path , opt.project_name )
	PCA_heatmap(df_DEGs,opt.project_name+'_DEGs', replicates_orders_list, opt.cuffdiff_path)
	Pearson_R(opt.cuffdiff_path, df_DEGs)
	PCA_heatmap(df_ALL_Genes, opt.project_name+'_All_Genes', replicates_orders_list,opt.cuffdiff_path)
	#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
