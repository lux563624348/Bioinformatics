"""
==================================
Template of the basics of violin plots
==================================

Violin plots are similar to histograms and box plots in that they show
an abstract representation of the probability distribution of the
sample. Rather than showing counts of data points that fall into bins
or order statistics, violin plots use kernel density estimation (KDE) to
compute an empirical distribution of the sample. That computation
is controlled by several parameters. This example demonstrates how to
modify the number of points at which the KDE is evaluated (``points``)
and how to modify the band-width of the KDE (``bw_method``).

For more information on violin plots and KDE, the scikit-learn docs
have a great section: http://scikit-learn.org/stable/modules/density.html
"""

import sys
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
import pandas as pd

def violin_plot(aaaa):
	
	File_Path=aaaa
	print File_Path
	
	df1=pd.read_excel(File_Path, sheetname=0, usecols=[0], convert_float=True, header=1)
	df2=pd.read_excel(File_Path, sheetname=0, usecols=[1], convert_float=True, header=1)
	
	df2=df2.dropna(axis=0, how='any')
	
	### generate two figure for comparing.
	fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

	df1.iloc[:,:]=np.log10(df1.iloc[:,:])
	df2.iloc[:,:]=np.log10(df2.iloc[:,:])

	# plot violin plot
	axes[0].violinplot(df1.iloc[:,0],
					   showmeans=False,
					   showmedians=True)
	axes[0].set_title('Solo Peaks Violin Plot')


	axes[1].violinplot(df2.iloc[:,0],
					   showmeans=False,
					   showmedians=True)
	axes[1].set_title('Common Peaks Violin Plot')
	# adding horizontal grid lines
	
	for ax in axes:
		ax.yaxis.grid(True)
		#ax.set_xticks([y + 1 for y in range(max(df1))])
		ax.set_xlabel('Frequency')
		ax.set_ylabel('Log10(Log10 pvalues)')


	plt.savefig('violin.png')
	plt.show()
	
	return 0


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="input_data_path", metavar="<file>", help="input_data")
	(opt, args)=parser.parse_args(argv)
	
	violin_plot(opt.input_data_path)
	if len(argv) < 3:
		parser.print_help()
	sys.exit(1)
	
	
	return 0
	
if __name__ == "__main__":
	main(sys.argv)

