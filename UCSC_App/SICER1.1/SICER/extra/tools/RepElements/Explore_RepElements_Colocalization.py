# Author :  Weiqun Peng


"""
Study the combinatorial structure of the RepElements
"""

import re, os, shutil, time, sys, operator
from math import *
from string import *
from optparse import OptionParser
import copy
from operator import itemgetter
try:
   import cPickle as pickle
except:
   import pickle

sys.path.append('/home/wpeng/data/SICER1.1/SICER/lib')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra')
sys.path.append('/home/wpeng/data/SICER1.1/SICER/extra/tools/RepElements')
sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')

import GenomeData 
import RepElements


plus = re.compile('\+');
minus = re.compile('\-');

RepElementsError = "Error in RepElement class";

def main(argv):

	parser = OptionParser()
	parser.add_option("-i", "--picklefile", action="store", type="string", dest="selected", metavar="<file>",
					help="pickle file for selected members in the format of {reClass:{reFamily:reName:[id]}}} ")    
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-n", "--name_for_all_pickle_files", action="store", type="string", dest="summary_name", help="common name of all pickle files, one pickle for one reName", metavar="<str>")
	
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	#load the selected ids
	print "Loading the pickle file for selected members in the format of {reClass:{reFamily:reName:[id]}}}"
	selected_members = pickle.load(open(opt.selected, 'rb'))
	
	
	
	# output the pkl file store the rep elements organized in a dictionary: {id: rep_element class instance}
	for reClass in mytree.keys():
		for reFamily in mytree[reClass].keys():
			for reName in mytree[reClass][reFamily]:
				re_file_name = "_".join([reClass, reFamily, reName]) + ".txt"
				known_repelements = RepElements.KnownRepElements.initiate_from_file(chroms, re_file_name)
				re_file_name = "_".join([reClass, reFamily, reName]) + ".pkl"
				known_repelements.output_pickle(re_file_name)
	

if __name__ == "__main__":
	main(sys.argv)
