#!/usr/bin/env/python
"""
Plot multiple curves from multiple files: 

input: 
1) a file that contains the list of file names for data used for plotting
	file_name   data_lable
	file_name	data_lable
	.....
2) the index of column used for plotting. The zeroth's column is for x values.

"""

import re, os, sys, pylab,shutil
from optparse import OptionParser

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/hg19/Annotation')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')


comment = re.compile("#|track")

def read_parameters(myfile):
	"""
	load the list of 
	[data file names]
	[data names: legends]
	"""
	file_names =[]
	legends = []
	f = open(myfile, 'r')
	for line in f:
		if not comment.match(line):
			line = line.strip()
			line = line.split()
			file_names.append(line[0])
			legends.append(line[1])
	return file_names, legends

def read_data(data_file, column_index=1):
	"""
	data file might have multiple columns, 
	the first column denotes the x values
	the second, third .. columns denotes the y1, y2, .. values
	
	return [x], [y]
	"""
	x = []
	y = []
	f = open(data_file,'r')
	for line in f:
		if not comment.match(line):
			line = line.strip()
			line = line.split()
			x.append(float(line[0]))
			y.append(float(line[column_index]))
	
	return x, y

 
def main(argv):
	parser=OptionParser()
	parser.add_option("-d","--file_for_data",action="store",type="string",dest="file_for_data",help="file that collects the name of data files to be plotted")
	parser.add_option("-c","--column_index", action="store",type="int",dest="column_index",help="column_index for y values")
	parser.add_option("-o","--outfile",action="store",type="string",dest="outfile",help="output file name ")
	(opt,args)=parser.parse_args()
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	#file_names: [file_name]
	#lengends: [legend str]
	file_names, legends = read_parameters(opt.file_for_data)
	print file_names
	print legends
	number_of_curves = len(file_names)
	
	pylab.figure(0)
	color = ['b', 'g', 'r', 'c', 'm', 'k']
	
	for i in xrange(number_of_curves):
		x, y = read_data(file_names[i], opt.column_index)
		pylab.plot(x, y, color[i], linewidth=3.0, label=legends[i])
	#pylab.xlabel(x_label,fontsize=30)
	pylab.ylabel('Normalized Read Counts' , fontsize=30)
	#pylab.title(title, fontsize=14)
	pylab.legend()	
	#pylab.show()
	pylab.savefig(opt.outfile + ".eps", format='eps')
	pylab.savefig(opt.outfile + ".png", format='png')
		
if __name__ == "__main__":
	main(sys.argv)

