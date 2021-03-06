#!/usr/bin/env python

import sys
import csv
from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqRecord import SeqRecord
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator

try:
   import cPickle as pickle
except:
   import pickle

sys.path.append('/home/data/hg19/JunZhu/modules/Resources/scripts')
sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RepElements')

import GenomeData 
import RepElements
import AnalyzeRNASeq

comment = re.compile ("#|track")

def get_sequences(summaries, genome, idset=None):
	"""
	summaries:[summary]
	
	summary pickle data structure 
	{id:
		{'annotation':RepElement class instance
		'index12_rc':
		'index12_rpkm':
		'index3_rc':
		'index3_rpkm':
		'index6_rc':
		'index6_rpkm':
		'index8_rc':
		'index8_rpkm':
		'index9_rc':
		'index9_rpkm':
		}
	}
	
	assembled REs: 
	{chrom: 
		{(region_start, region_end):
			{"elements":[ids]; 
			"age":  
			target_name + "_rc":
			control_name + "_rc":
			"strand": "+";
			"num_boundary_elements": value; 
			"5_boundary_elements":[id]; 
			"5_boundary_elements_age":
			"3_boundary_elements":[id];
			"3_boundary_elements_age":
			"I_boundary_elements":[id]; 
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	
	
	return {id:SeqRecord object}
	
	"""
	#get subset
	mysummaries = AnalyzeRNASeq.get_subset(summaries, idset)
	
	#Separate by chrom
	mysummaries_by_chrom = [] #[{chrom:{id:{feature_name:value}}}]
	for mysummary in mysummaries:
		mysummaries_by_chrom.append(AnalyzeRNASeq.separate_by_chrom(mysummary))
		
	myseqs = {}
	for chrom_seq_record in SeqIO.parse(open(genome, "rU"), 'fasta'): #iterate through chrom first to save memory
		chrom = chrom_seq_record.id
		for mysummary_by_chrom in mysummaries_by_chrom:#{chrom:{id:{feature_name:value}}}
			if chrom in mysummary_by_chrom.keys():
				mysummary_this_chrom = mysummary_by_chrom[chrom]
				#print chrom, len(mysummary_this_chrom.keys())
				for myid in mysummary_this_chrom.keys():
					myre = mysummary_this_chrom[myid]['annotation']
					mystrand = myre.strand
					mystart = myre.genoStart
					myend = myre.genoEnd
					re_seq = chrom_seq_record[mystart:myend] #SeqRecord object: id, seq, etc
					re_seq.id = myid
					if mystrand == "-":
						re_seq.seq = (chrom_seq_record[mystart:myend].seq).reverse_complement() #Seq object
					myseqs[myid] = re_seq
	#print "The number of elements is %d" %len(myseqs)
	return myseqs
	
def get_LTR_clusters(assembled_REs):
	"""
	Collect the LTR information from the assembled_REs
	the strand information are not from the LTR themselves, but from the int parts
	return {(id):strand}
	"""
	RE_LTR_clusters = {}
	for chrom in assembled_REs:
		for myregion in assembled_REs[chrom]:
			strand = assembled_REs[chrom][myregion]["strand"]
			five_LTR_elements = assembled_REs[chrom][myregion]["5_boundary_elements"] #[id]
			if len(five_LTR_elements) > 0 :
				assert strand == AnalyzeRNASeq.find_strand_from_id(five_LTR_elements[0])
				RE_LTR_clusters[tuple(five_LTR_elements)] = strand
			three_LTR_elements = assembled_REs[chrom][myregion]["3_boundary_elements"] #[id]
			if len(three_LTR_elements) > 0 :
				assert strand == AnalyzeRNASeq.find_strand_from_id(three_LTR_elements[0])
				RE_LTR_clusters[tuple(three_LTR_elements)] = strand
	return RE_LTR_clusters

def get_coordinates(summaries, idset):
	"""
	summaries:[summary]
	
	summary pickle data structure 
	{id:
		{'annotation':RepElement class instance
		'index12_rc':
		'index12_rpkm':
		'index3_rc':
		'index3_rpkm':
		'index6_rc':
		'index6_rpkm':
		'index8_rc':
		'index8_rpkm':
		'index9_rc':
		'index9_rpkm':
		}
	}
	
	return: {id:(strand, start, end)}
	
	"""
	
	coordinates = { }
	#get subset
	mysummaries = AnalyzeRNASeq.get_subset(summaries, idset)
	
	for summary in mysummaries:
		for myid in summary.keys():
			myre = summary[myid]['annotation']
			mystrand = myre.strand
			mystart = myre.genoStart
			myend = myre.genoEnd
			coordinates[myid] = (mystrand, mystart, myend)
	return coordinates

def get_coordinates_of_assembled_REs(summaries, assembled_REs):
	"""
	Collect the coordinate information from the assembled_REs: from start of 5LTR to end of 3LTR
	the strand information are not from the LTR themselves, but from the int parts
	return {chrom :[(ids, strand, (start, end), #LTR]}
	"""
	assembled_RE_coordinates = {}
	
	for chrom in assembled_REs:
		assembled_RE_coordinates[chrom] = []
		for myregion in assembled_REs[chrom]:
			num_LTR = 0
			int_ids = assembled_REs[chrom][myregion]["elements"]
			int_start = myregion[0]
			int_end = myregion[1]
			assert (int_start <= int_end)
			strand = assembled_REs[chrom][myregion]["strand"]
			
			five_LTR_elements = assembled_REs[chrom][myregion]["5_boundary_elements"] #[id]
			#print "five_LTR_elements", five_LTR_elements
			five_LTR_start = int_start
			five_LTR_end = int_end
			if len(five_LTR_elements) > 0 :
				num_LTR += 1
				
				#{id:(strand, start, end)}
				coordinates = get_coordinates(summaries, five_LTR_elements)
				#print coordinates
				LTR_starts = [item[1] for item in coordinates.values()]
				LTR_ends = [item[2] for item in coordinates.values()]
				five_LTR_start = min(LTR_starts)
				five_LTR_end = max(LTR_ends)
				assert strand == AnalyzeRNASeq.find_strand_from_id(five_LTR_elements[0])
				
			
			three_LTR_elements = assembled_REs[chrom][myregion]["3_boundary_elements"] #[id]
			#print "three_LTR_elements", three_LTR_elements
			three_LTR_end = int_end
			three_LTR_start = int_start
			if len(three_LTR_elements) > 0 :
				num_LTR += 1
				
				coordinates = get_coordinates(summaries, three_LTR_elements)
				LTR_starts = [item[1] for item in coordinates.values()]
				LTR_ends = [item[2] for item in coordinates.values()]
				three_LTR_start = min(LTR_starts)
				three_LTR_end = max(LTR_ends)
				assert strand == AnalyzeRNASeq.find_strand_from_id(three_LTR_elements[0])
			
			ids = tuple(five_LTR_elements + int_ids + three_LTR_elements)
			if strand == "+":
				assert (five_LTR_start <= int_start)
				assert (three_LTR_end >= int_end)
				assembled_RE_coordinates[chrom].append((ids, strand, (five_LTR_start, three_LTR_end), num_LTR))
			elif strand == "-":
				assert (three_LTR_start <= int_start)
				assert (five_LTR_end >= int_end)
				assembled_RE_coordinates[chrom].append((ids, strand, (three_LTR_start, five_LTR_end), num_LTR))
			
	return assembled_RE_coordinates
	
def get_sequences_of_assembled_REs(summaries, assembled_REs, genome, num_LTR):
	"""
	Get the sequences of the assembled_REs, including "5_boundary_elements", int_elements,  "3_boundary_elements"
	return {id:sequence}
	"""
	# {chrom :[(ids, strand, (start, end), #LTR]}
	assembled_RE_coordinates = get_coordinates_of_assembled_REs(summaries, assembled_REs)
	
	myseqs = {}
	for chrom_seq_record in SeqIO.parse(open(genome, "rU"), 'fasta'): #iterate through chrom first to save memory
		chrom = chrom_seq_record.id
		if chrom in assembled_RE_coordinates.keys():
			for item in assembled_RE_coordinates[chrom]:
				
				my_num_LTR = item[3]
				if my_num_LTR == num_LTR:
					ids = item[0]
					myid = ids[0] # use the first id to represent the whole thing
					mystrand = item[1]
					
					mystart = item[2][0]
					myend = item[2][1]
					
					
					re_seq = chrom_seq_record[mystart:myend] #SeqRecord object: id, seq, etc
					re_seq.id = ids[0]
					if mystrand == "-":
						re_seq.seq = (chrom_seq_record[mystart:myend].seq).reverse_complement() #Seq object
					myseqs[myid] = re_seq
	return myseqs
	
	
def assemble_sequences(myseqs, RE_LTR_clusters):
	"""
	The functional unit for a RE is not a single RE, rather it is a cluster of REs.
	myseqs: [SeqRecord object]
	RE_LTR_clusters: derived from asembled_RE pickle see AnalyzeRNASeq.py: {[id]:strand}
	return {id:SeqRecord object}
	"""
	my_assembled_seqs = {}
	for cluster in RE_LTR_clusters.keys(): #(id)
		strand = RE_LTR_clusters[cluster]
		myids = list(set(myseqs.keys()) & set(cluster))
		if len(myids) > 0:
			my_assembled_seq = Seq("")
			my_assembled_seq_id = ""
			if strand == "+":
				myids.sort()
				for myid in myids:
					my_assembled_seq += myseqs[myid].seq
					my_assembled_seq_id += myid + ","
			elif strand == "-": 
				myids.sort(reverse=True)
				for myid in myids:
					my_assembled_seq += myseqs[myid].seq #each of them records the reverse complement
					my_assembled_seq_id += myid + ","
			
			#print "This is", my_assembled_seq_id
			my_assembled_seqs[my_assembled_seq_id] = SeqRecord(my_assembled_seq, id=my_assembled_seq_id, description="")
	return my_assembled_seqs
	
def main(argv):

	parser=OptionParser()
	parser.add_option("-i", "--pickle_file_name_for_RE", action="store", type="string",
						dest="summary_pickle", help="summary pickle for a particular RE", metavar="<file>")		
	parser.add_option("-g", "--genome", action="store", type="string",
						dest="genome", help="genome", metavar="<file>")
	parser.add_option("-t", "--id_subset", action="store", type="string",
						dest="id_subset_file", help="id subset file", default="ALL", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
						dest="output_file_name", help="output sequence file name", metavar="<file>")					
	(opt,args)=parser.parse_args()
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	#load the summary file
	print "Loading ", opt.summary_pickle
	summary = pickle.load(open(opt.summary_pickle, 'rb'))
	
	if opt.id_subset_file == "ALL":
		id_subset = None
	else:
		id_subset = []
		inf = open(opt.id_subset_file, "r")
		for line in inf:
			if not re.match("#", line):
				line = line.strip();
				sline = line.split("\t");
				id_subset.append(sline[0])
	
	names =  AnalyzeRNASeq.find_family_name([summary])
	print "The number of elements in %s is %d" %(names[0], len(summary.keys()))
	myseqs = get_sequences([summary], opt.genome, id_subset) #{id:SeqRecord object}}
	print "The number of sequences is %d" %len(myseqs.keys())
	myseqs_list = myseqs.values()
	output_handle = open(opt.output_file_name, "w")
	SeqIO.write(myseqs_list, output_handle, "fasta")
	output_handle.close()
	
	
	
	
	#Get the sequence of the assembled REs, to 
	"""
	assembled REs: 
	{chrom: 
		{(region_start, region_end):
			{"elements":[ids]; 
			"age":  
			target_name + "_rc":
			control_name + "_rc":
			"strand": "+";
			"num_boundary_elements": value; 
			"5_boundary_elements":[id]; 
			"5_boundary_elements_age":
			"3_boundary_elements":[id];
			"3_boundary_elements_age":
			"I_boundary_elements":[id]; 
			"expression_fc_"+target_name+"_vs_" + control_name: max_value
			}
		}
	}
	"""
	
	current_dir = os.getcwd()
	path = "/home/data/mm9/Lin/processed/RepElements/brokendown/summary"
	os.chdir(path)
	print "Loading the Assembled_REs.pkl"
	assembled_REs = pickle.load(open("Assembled_REs.pkl", 'rb'))
	num_re = 0
	for chrom in assembled_REs.keys():
		num_re += len(assembled_REs[chrom])
	print "There are %d elements in Assembled_REs.pkl" %num_re
	os.chdir(current_dir)
	
	
	##Get the LTR clustering information from Assembled_REs.pkl, the strand of the LTR is determined from the int.
	#RE_LTR_clusters = get_LTR_clusters(assembled_REs) #{(id):strand}
	#print "The number of LTR clusters in Assembled_REs.pkl is ", len(RE_LTR_clusters)
	#my_assembled_seqs =  assemble_sequences(myseqs, RE_LTR_clusters)
	#print "The number of assembled elements is %d" %len(my_assembled_seqs.keys())
	
	#my_assembled_seqs_list = my_assembled_seqs.values()
	#output_handle = open(opt.output_file_name + "_assembled", "w")
	#SeqIO.write(my_assembled_seqs_list, output_handle, "fasta")
	#output_handle.close()
	
	
	#Get the full assembled_RE sequences
	current_dir = os.getcwd()
	path = "/home/data/mm9/Lin/processed/RepElements/brokendown/summary"
	os.chdir(path)
	summaries_names = ["summary_on_LTR_ERV1_RLTR4_MM-int.pkl", "summary_on_LTR_ERV1_MuLV-int.pkl", "summary_on_LTR_ERV1_RLTR4_Mm.pkl"]
	summaries = []
	for summary_file_name in summaries_names:
		inf = open(summary_file_name, 'rb')
		summaries.append (pickle.load(inf))
		inf.close()
	os.chdir(current_dir)
	
	# Get the 2-LTRs
	#{id:sequence}
	for num_LTR in [0,1,2]:
		my_assembled_seqs = get_sequences_of_assembled_REs(summaries, assembled_REs, opt.genome, num_LTR)
		my_assembled_seqs_list = my_assembled_seqs.values()
		print "The number of %d LTR elements is %d." %(num_LTR, len(my_assembled_seqs_list))
		output_handle = open("assembled_" + str(num_LTR) + "LTR_RLTR4_sequences.fa", "w")
		SeqIO.write(my_assembled_seqs_list, output_handle, "fasta")
		output_handle.close()

if __name__ == "__main__":
    main(sys.argv)

