## INTRODUCTION
########################################################################
## 11/1/2018
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Ver.beta
########################################################################

import numpy as np
import pandas as pd
from optparse import OptionParser
import sys

### Filtering 
def Return_filtered_peaks_bed_format(path_folder, name, fold_change, pvalue, qvalue):
	df = pd.read_csv(path_folder+'/'+name, sep='\t', header=0, comment='#')
	number_peaks = len(df)
	df = df[ (df['-log10(pvalue)']>-np.log10(pvalue)) & (df['-log10(qvalue)']>-np.log10(qvalue)) & (df['fold_enrichment']>fold_change )].loc[:,['chr','start','end','name']]
	df['name'] = df['name'].str[-10:]
	df = df.rename(columns={'chr':'#chr'})
	df.to_csv(path_folder + '_'+str(number_peaks)+ '_filtered_peask_'+ str(len(df))+'.bed', sep='\t', index=False, header=True)
	return df


def Return_peaks_bed_format(path_folder, name):
	df = pd.read_csv(path_folder+'/'+name, sep='\t', header=0, comment='#')
	number_peaks = len(df)
	df['name'] = df['name'].str[-10:]
	df = df.loc[:,['#chr','start','end','name']].to_csv(path_folder+'/'+name[:-4]+ '_'+str(number_peaks)+'.bed', sep='\t', index=True, header=True)
	return df


def main(argv):
	desc="""This is a template for the Macs2 peaks output filtering."""
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--input_folder", action="store", type="string",
			dest="input_folder", help="Input Folder")
	parser.add_option("-n", "--name", action="store", type="string",
			dest="input_name", help="Input Name", metavar="<file>")
	parser.add_option("-f", "--fold_change", action="store", type="float", dest="fc",
			help="fold change shall be large than this value")
	parser.add_option("-p", "--p_value", action="store", type="float", dest="p_value",
			help="p value shall be less than this value")
	
	parser.add_option("-q", "--q_value", action="store", type="float", dest="q_value",
			help="q value shall be less than this value")
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 5:
		parser.print_help()
		sys.exit(1)
	
	print " "
	print "Here is the Summary of your input."
	print "Name of Input File: %s" % opt.input_name
	print "Fold Change threshold: " + str(opt.fc)
	print "p value threshold " + str(opt.p_value)
	print "q value threshold " + str(opt.q_value)
	print " "
	
	Return_filtered_peaks_bed_format(opt.input_folder, opt.input_name, opt.fc, opt.p_value, opt.q_value)
	
	return 0
	
if __name__ == "__main__":
	main(sys.argv)
