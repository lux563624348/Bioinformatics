#!/usr/bin/env/python
import re, os, sys, pylab,shutil
from optparse import OptionParser

def read_data_from_file(datafile,x=[],readp=[],readm=[]):
    for line in file(datafile):
        a,b,c=map(float,line.split())
        x.append(a)
        readp.append(b)
        readm.append(c)

def plot_data(x,read,read_b,fignum,range,title,xtitle,legend,outgraphname):
    pylab.figure(fignum)
    pylab.plot(x,read,"b",label=legend[0])
    pylab.plot(x,read_b,'r',label=legend[1])
    pylab.xlabel(xtitle,fontsize=14)
    pylab.xlim((-1*range),range)
    pylab.ylabel('Normalized Counts' , fontsize=14)
    pylab.title(title)
    pylab.legend()
    pylab.show()
    pylab.savefig(outgraphname)
    
    
def main(argv):
    parser=OptionParser()
    parser.add_option("-a","--input_graph_file_a",action="append",type="string",dest="graphfile",help="input graph file a")
    parser.add_option("-b","--input_graph_file_b",action="append",type="string",dest="graphfile",help="input graph file b")
    parser.add_option("-c","--legend_a",action="append",type="string",dest="legend",help="legend string for file a")
    parser.add_option("-d","--legend_b",action="append",type="string",dest="legend",help="legend string for file b")
    parser.add_option("-t","--title",action="store",type="string",dest="title",help="plot title")
    parser.add_option("-x","--xtitle",action="store",type="string",dest="xtitle",help="plot x axis title")
    parser.add_option("-r","--range",action="store",type="float",dest="range",help="symmetric range for graph")
    parser.add_option("-o","--outfile",action="store",type="string",dest="outfile",help="legend string for graph file a")

    (opt,args)=parser.parse_args()
    if len(argv) < 16:
        	parser.print_help()
        	sys.exit(1)
	
    x_a=[]
    readp_a=[]
    readm_a=[]
    read_data_from_file(opt.graphfile[0],x_a,readp_a,readm_a)

    x_b=[]
    readp_b=[]
    readm_b=[]
    read_data_from_file(opt.graphfile[1],x_b,readp_b,readm_b)

    plot_data(x_a,readp_a,readp_b,0,opt.range,opt.title+' (+ reads)',opt.xtitle,opt.legend,opt.outfile+'_preads')
    plot_data(x_a,readm_a,readm_b,1,opt.range,opt.title+' (- reads)',opt.xtitle,opt.legend,opt.outfile+'_nreads')
   

if __name__ == "__main__":
    main(sys.argv)

