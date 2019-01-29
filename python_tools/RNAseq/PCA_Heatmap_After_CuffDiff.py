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

import numpy as np
import pandas as pd
import os
import sys
from optparse import OptionParser
### Parameters for DEGs:

FPKM_threshold=1.0
####################################################################################


## FUNCTIONS
### FUNCTION
# For any input of /gene_exp.diff, return its All genes with
### gene_id, cond1, cond2, log2(fold_change), p_value, plus a number of order.
def generate_All_Genes(Input_Path, number):
#### READ FILE
    df = pd.read_csv(Input_Path+'/gene_exp.diff', sep='\t', header=0, usecols={'test_id',\
    'status','sample_1','sample_2','value_1','value_2','log2(fold_change)','p_value','q_value'})
#### Rename columns
    df=df.rename(columns={'test_id':'gene_id', 'value_1': df['sample_1'].unique()[0], 'value_2': df['sample_2'].unique()[0]})
#### Output   
    return df.loc[:,['gene_id',df['sample_1'].unique()[0],df['sample_2'].unique()[0],'log2(fold_change)',
                     'p_value','q_value', 'End_'+str(number), '||']].fillna('')
####################################################################################
def generate_Upregulated_Genes(Input_Path):
#### READ FILE
    df = pd.read_csv(Input_Path+'/gene_exp.diff', sep='\t', header=0, usecols={'test_id',\
    'status','sample_1','sample_2', 'value_1','value_2','log2(fold_change)','p_value','q_value'})
#### Filter
    df=df[(df['status']=='OK') & (df['q_value']<=q_value_less) & (df['value_2']>=FPKM_threshold) & (df['log2(fold_change)'] >= np.log2(FC_UP))]
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
#### Functions to read /genes.read_group_tracking and output FPKM
def generate_genes_FPKM_df(Input_Path):
    df = pd.read_csv(Input_Path+'/genes.read_group_tracking', sep='\t', header=0)
    FPKM_df = None
    for cond in df['condition'].unique():
        for replica in df[df['condition'] == cond]['replicate'].unique():
            #print (cond, replica)
            temp_df = (df[ (df['condition'] == cond) & (df['replicate'] == replica)].loc[:,['tracking_id', 'FPKM']])
            cond_name= cond+'_'+str(replica+1)
            temp_df.rename(columns={'tracking_id': 'gene_id', 'FPKM': cond_name }, inplace=True)
            if FPKM_df is None:
                FPKM_df
                FPKM_df = temp_df  
            else:
                FPKM_df = FPKM_df.merge(temp_df, on='gene_id', how='outer', suffixes=('','_'))
    return FPKM_df
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
####################################################################################
def Summary_Expression(Cuffdiff_Path, output_name):
    import datetime
    PATH_FOLDER=Cuffdiff_Path+ '/CuffDiff_Results'
    OUT_FOLDER=Cuffdiff_Path+'/CuffDiff_Results_Summary'

    DIR_CHECK_CREATE(PATH_FOLDER)
    DIR_CHECK_CREATE(OUT_FOLDER)

    INPUT_LIST=[f for f in os.listdir(PATH_FOLDER) if not f.startswith('.')]
    writer = pd.ExcelWriter(OUT_FOLDER+'/'+output_name+'_CuffDiff_Summary_'+str(datetime.date.today())+ '.xlsx', engine='xlsxwriter')

    i=0
    for input_name in INPUT_LIST[:]:
        INPUT_PATH = PATH_FOLDER+'/'+input_name
        print INPUT_PATH
        if (i==0):
            df_all=generate_genes_FPKM_df(INPUT_PATH)
            df_all[input_name]=''
            df_all = df_all.merge(generate_All_Genes(INPUT_PATH,i+1), on='gene_id', how='inner', suffixes=('','_')) 
            i+=1
            continue
        df_all = df_all.merge(generate_genes_FPKM_df(INPUT_PATH), on='gene_id', how='inner', suffixes=('','_'))
        df_all[input_name]=''
        df_all = df_all.merge(generate_All_Genes(INPUT_PATH,i+1), on='gene_id', how='inner', suffixes=('','_')) 

        i+=1
    df_all.to_excel(writer, sheet_name='All_Genes', index=None)

    ###### Generating a DEGs list in the last.
    DEGs_List=pd.DataFrame([],columns=['gene_id'])
    
    for input_name in INPUT_LIST:
        INPUT_PATH = PATH_FOLDER+'/'+input_name
        print ('Library:' + input_name)
        df_up = generate_genes_FPKM_df(INPUT_PATH)
        df_up[input_name]=''
        df_up = df_up.merge(generate_Upregulated_Genes(INPUT_PATH), on='gene_id', how='inner', suffixes=('','_'))
        df_up.to_excel( writer, sheet_name='up_'+input_name, index=None)

        DEGs_List=pd.concat([DEGs_List,df_up[['gene_id']]])
        print ('# of Up:' )
        print(df_up.shape)

        df_down = generate_genes_FPKM_df(INPUT_PATH)
        df_down[input_name]=''
        df_down = df_down.merge(generate_Downregulated_Genes(INPUT_PATH), on='gene_id', how='inner', suffixes=('','_')) 
        df_down.to_excel( writer, sheet_name='down_'+input_name,index=None)
        DEGs_List=pd.concat([DEGs_List,df_down[['gene_id']]])
        print ('# of Down')
        print(df_down.shape)
        print ('')

    print ("Total Number of DEGs is:" + str(len(DEGs_List['gene_id'].unique())))

    DEGs_List = pd.DataFrame( list(DEGs_List['gene_id'].unique()), columns=['gene_id'])
    DEGs_List.to_excel( writer, sheet_name='Union_DGEs_List',index=None)
    writer.save()
    ### Output A Excel Summary Completed!
    All_Gene_List = pd.DataFrame( list(df_all['gene_id'].unique()), columns=['gene_id'])
    
    df_DEGs=DEGs_List
    df_All_Genes=All_Gene_List
    for input_name in INPUT_LIST[:]:
        INPUT_PATH = PATH_FOLDER+'/'+input_name
        print ('Library:' + input_name)
        df_DEGs = df_DEGs.merge(generate_genes_FPKM_df(INPUT_PATH), on='gene_id',how='inner')
        df_All_Genes = df_All_Genes.merge(generate_genes_FPKM_df(INPUT_PATH), on='gene_id',how='inner')
    print ('# of DEGs:' )
    print(df_DEGs.shape)
    print ('# of All Genes')
    print(df_All_Genes.shape)
  
    
    return df_DEGs.set_index('gene_id'),df_All_Genes.set_index('gene_id')
####################################################################################
def PCA_heatmap(df, title_name, lib_order, Cuffdiff_Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    OUT_FOLDER=Cuffdiff_Path+'/CuffDiff_Results_Summary'
    sns.set_style("white") 
    ### Drop all zeros rows for plot
    Heatmap_Df = df.loc[(df!=0).any(axis=1)]
    print ('# of Genes (Excluding genes that RPKM of all conditions equal zero.')
    print(Heatmap_Df.shape)
    Heatmap_Name='Heatmap_'+title_name +'_'+ str(len(Heatmap_Df.index))
    #### Very Important
    #Heatmap_Df=FPKM_df.fillna(0)
    fig_1 = sns.clustermap(Heatmap_Df,  yticklabels=False, z_score=0, col_cluster=True, cmap='RdBu_r' )
    fig_1.fig.suptitle(Heatmap_Name)
    fig_1.savefig(OUT_FOLDER+'/'+Heatmap_Name+'.png')
    print (OUT_FOLDER+'/'+Heatmap_Name+'.png  generated!' )
    
    ###### PCA
    # PCA
    df_2 = df #.reset_index()

    #### Remove unicode of list in python 
    df2_gene_id = [x.encode('ascii', 'ignore') for x in df_2.index]


    from sklearn.decomposition import PCA
    from sklearn import datasets
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.preprocessing import StandardScaler

    #### In here we set targets as the name of columns, which means our purpose is to compare \
    #### the relationship between different columns.
    targets= [x.encode('ascii', 'ignore') for x in df_2.columns]

    ################################################################################

    df_2_T = df_2.transpose()

    # Separating out the features
    x = df_2_T.loc[ :, df2_gene_id ].values


    # Separating out the target
    y = df_2_T.loc[targets,:].values


    # Standardizing the featuresbio
    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ['principal component 1', 'principal component 2'])


    ########### THIS is fucking killing me. reindex df with df2, just using following command.
    principalDf.index = df_2_T.index
    
    order=lib_order # Each lib has at most 8 replicates so far.
    ### PLOT
    ### Some Parameters
    matplotlib.rcParams['font.sans-serif'] = ['Arial'] #, 'sans-serif']
    fig = plt.figure(figsize = (10,8))
    ax = fig.add_subplot(1,1,1) 
    s_size=300
    alpha_value=1
    #### In here we set targets as the name of columns, which means our purpose is to compare \
    #### the relationship between different columns.
    targets= [x.encode('ascii', 'ignore') for x in df_2.columns]
    filled_markers = ( 'o','^','*','X', 'D', 's', 'd', 'P' ) 
    if len(targets) == sum(order):
        colors=[]
        markers=[]
        k=0
        for i in order[:]: 
            for j in range(i):
                colors.append("C"+ str(k))
                markers.append(filled_markers[j])
            k+=1
    ### Associated each replicate with color and marker

    for target, color, mark_er in zip( targets,colors,markers):
        indicesToKeep = principalDf.index == target
        ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1']
                   , principalDf.loc[indicesToKeep, 'principal component 2']
                   , s = s_size, c=color, alpha=alpha_value, marker=mark_er, label=target)

    legend_object = ax.legend(loc="upper right", bbox_to_anchor=(0.4, 1, 1,0),edgecolor='w',
              borderaxespad=0,fancybox=True, shadow=False,  fontsize=16)    

    # change the font colors to match the line colors:
    i=0
    for text in legend_object.get_texts():
        text.set_color(colors[i])
        i+=1
    #change the font colors to match the line colors:    

    ax2 = ax.twiny()
    ax3 = ax.twinx()
    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax3.yaxis.set_ticks_position("left")
    ax3.yaxis.set_label_position("left")

    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.05))
    ax2.spines["bottom"].set_linewidth(2)
    ax3.spines["left"].set_position(("axes", -0.05))
    ax3.spines["left"].set_linewidth(2)


    #ax2.set_frame_on(True)
    #ax3.set_frame_on(True)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])

    ax2.set_xlim(ax.get_xlim())
    ax3.set_ylim(ax.get_ylim())
    ax2.tick_params(axis='x',which='major', direction='out', length=6, labelsize=18)
    ax2.tick_params(axis='x',which='minor', direction='out', length=6, labelsize=18)
    ax3.tick_params(axis='y',which='major', direction='out', length=6, labelsize=18)
    ax3.tick_params(axis='y',which='minor', direction='out', length=6, labelsize=18)

    ax2.set_xlabel('Principal Component 1', fontsize = 22, fontname='Arial')
    ax3.set_ylabel('Principal Component 2', fontsize = 22, fontname="Arial")

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    
    fig.tight_layout()
    fig.savefig( OUT_FOLDER+'/'+'PCA On '+title_name+'.png')
    print (OUT_FOLDER+'/'+'PCA On '+title_name+'.png  generated!' )

    return None
####################################################################################
def Pearson_R(Cuffdiff_Path,df):
    from scipy.stats import pearsonr
    OUT_FOLDER=Cuffdiff_Path+'/CuffDiff_Results_Summary'
    R_output_matrix=pd.DataFrame(columns=df.columns, index=df.columns)
    Pvalue_output_matrix=pd.DataFrame(columns=df.columns, index=df.columns)
    for i in range(len(df.columns)):
        for j in range(len(df.columns)):
            R_output_matrix.iloc[i,j]= (round(pearsonr( (df.iloc[:,i]), (df.iloc[:,j]))[0], 3), round(pearsonr( (df.iloc[:,i]), (df.iloc[:,j]))[1], 3))
    R_output_matrix.to_csv(OUT_FOLDER+'/Pearson_Correlation_Summary.txt', sep='\t')
    print " "
    print ("Correlation Summary Of All Replicates")
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
