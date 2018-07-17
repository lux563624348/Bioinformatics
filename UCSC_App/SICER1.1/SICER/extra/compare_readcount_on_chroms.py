#!/usr/bin/env python

# Calculate read counts on chroms


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

def get_rc(rc_file):
	"""
	"# Name"  + "\t" + "RPM" + "\t"  + "RPKM" + '\n'
	rc_info: {name:[rc, rpm, rpkm]}
	
	"""
	rc_info = {}
	
	inf = open(rc_file,"r")
	for line in inf:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			chrom = sline[0]
			rc_info[chrom] = [float(sline[i]) for i in range(1,len(sline)) ] 
	inf.close()
	return rc_info
	
def get_total_rc(rc_info):
	for name in rc_info.keys():
		rc = rc_info[name][0]
		if rc > 5:
			rpm = rc_info[name][1] #rpm = rc/total_tc * 1000000
			total_rc = rc * 1000000.0/ rpm 
			return total_rc
	
def combine_and_calculate_fc(control_dic, target_dic, pseudo_count=1):
	"""
	control_dic: {name:[rc, rpm, rpkm]}
	target_dic: {name:[rc, rpm, rpkm]}
	"""
	combined_dic={}
	shared_keys = list(set(control_dic.keys()) & set (target_dic.keys()))
	control_dic_keys  = list(set(control_dic.keys()) - set(shared_keys))
	target_dic_keys  = list(set(target_dic.keys()) - set(shared_keys))
	
	total_rc_control = get_total_rc(control_dic)
	total_rc_target = get_total_rc(target_dic)
	
	for mykey in shared_keys:
		fc = ((target_dic[mykey][0] + pseudo_count)/total_rc_target)/((control_dic[mykey][0] + pseudo_count)/total_rc_control)
		combined_dic[mykey] = (target_dic[mykey][0], target_dic[mykey][2], control_dic[mykey][0], control_dic[mykey][2],  fc)
	for mykey in control_dic_keys:
		fc = ((0 + pseudo_count)/total_rc_target)/((control_dic[mykey][0] + pseudo_count)/total_rc_control)
		combined_dic[mykey] = (0, 0, control_dic[mykey][0], control_dic[mykey][2], fc)
	for mykey in target_dic_keys:
		fc = ((target_dic[mykey][0] + pseudo_count)/total_rc_target)/((0 + pseudo_count)/total_rc_control)
		combined_dic[mykey] = (target_dic[mykey][0], target_dic[mykey][2], 0, 0,  fc)
	return combined_dic
		
def output_dic(dic, filename):
	f = open(filename,"w")
	outline = "#Name" + "\t" + "Target_RC" + "\t"  + "Target_RPKM" + "\t" + "Control_RC" + "\t"  + "Control_RPKM" + "\t" + "Target_vs_control_fc" + '\n'
	f.write(outline)
	sorted_keys = sorted(dic.keys())
	for key in sorted_keys:
		mylist = dic[key]
		mylist_str = "\t".join([str(i) for i in mylist])
		outline = str(key) + "\t" + mylist_str + "\n"
		f.write(outline)
	f.close()
		
def main(argv):
	parser = OptionParser()
	
	parser.add_option("-i", "--TargetRcfile", action="store", type="string", dest="rc_file_target", metavar="<file>", help="target read count file ")
	parser.add_option("-j", "--ControlRcfile", action="store", type="string", dest="rc_file_control", metavar="<file>", help="control read count file ")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	target_dic = get_rc(opt.rc_file_target) #target_dic: {name:[rc, rpm, rpkm]}
	control_dic = get_rc(opt.rc_file_control) #control_dic: {name:[rc, rpm, rpkm]}
	
	pseudo_count = 1
	#
	combined_dic = combine_and_calculate_fc(control_dic, target_dic, pseudo_count)
	
	output_dic(combined_dic, opt.out_file)

if __name__ == "__main__":
	main(sys.argv)