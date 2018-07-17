#!/usr/bin/env/python
# Plot two curves from a single file: x from the first column, y1 from the second column, y2 the third column

import re, os, sys, pylab,shutil
from optparse import OptionParser

def read_data_from_file(datafile,x=[],readp=[],readm=[]):
	for line in file(datafile):
		a,b,c=map(float,line.split())
		x.append(a)
		readp.append(b)
		readm.append(c)

def plot_data(x,readp, readm, myrange, title, x_label, legend, outgraphname):
	pylab.figure(0)
	line1, line2 = pylab.plot(x,readp,"b-", x, readm, "r-")
	pylab.legend((line1, line2), (legend[0], legend[1]), 'upper right' )
	pylab.xlabel(x_label,fontsize=14)
	pylab.xlim((-1*myrange),myrange)
	pylab.ylabel('Normalized Counts' , fontsize=14)
	pylab.title(title, fontsize=10)
	#pylab.show()
	pylab.savefig(outgraphname, format='eps')
    
    
def main(argv):
	parser=OptionParser()
	parser.add_option("-a","--data_file ",action="store",type="string",dest="datafile",help="data file of format x, y1, y2")
	parser.add_option("-b","--legend1",action="store",type="string",dest="legend1",help="legend string for y1")
	parser.add_option("-c","--legend2",action="store",type="string",dest="legend2",help="legend string for y2")
	parser.add_option("-t","--title",action="store",type="string",dest="title",help="plot title")
	parser.add_option("-x","--x_label",action="store",type="string",dest="x_label",help="plot x axis label")
	parser.add_option("-r","--myrange",action="store",type="float",dest="myrange",help="symmetric range for graph")
	parser.add_option("-o","--outfile",action="store",type="string",dest="outfile",help="output file name ")

	(opt,args)=parser.parse_args()
	if len(argv) < 14:
		parser.print_help()
		sys.exit(1)
	
	x=[]
	readp=[]
	readm=[]
	
	#print opt.datafile
	# print opt.x_label
	#print opt.myrange
	
	read_data_from_file(opt.datafile, x, readp, readm)
	legends=[opt.legend1, opt.legend2]

	plot_data(x, readp, readm, opt.myrange, opt.title, opt.x_label, legends, opt.outfile)

if __name__ == "__main__":
	main(sys.argv)

