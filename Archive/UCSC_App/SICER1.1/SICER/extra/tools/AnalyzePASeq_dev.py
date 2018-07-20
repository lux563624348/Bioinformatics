#!/usr/bin/env python
# Build a data structure: gene_id: [list of PA peaks]
# Calculate the characteristics for each gene_id

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
from operator import itemgetter
try:
   import cPickle as pickle
except:
   import pickle
from operator import itemgetter

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')
sys.path.append('/home/data/SICER1.1/SICER/extra/tools/RNASeq')
sys.path.append('/home/data/hg19/Annotation')

import Utility_extended
import GenomeData 
import SeparateByChrom
import UCSC_revised
import Entrez

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")


def AssignPeaksToEntrez3UTRs(entrez_genes, peakfile, chroms, chrom_lengths, peak_threshold, downstream_extension):
	"""
	Returns {entrez_id:(gene, ThreeUTR_length, peaks_on_3UTR)}
	gene:gene = entrez_genes_by_chrom.entrez_genes[entrez_id] 
	ThreeUTR_length: longest 3UTR length; length includes the downstream extension
	peaks_on_3UTR:[(location, read_count)]
	"""
	
	peaks_on_entrez_3UTRs = {} #store the peaks for each 3UTR of the entrez cluster. {Entrez_ID: (gene, ThreeUTR_length, peaks_on_3UTR)}
	
	if Utility_extended.fileExists(peakfile):
		# Read the peaks, which is assumed to have the pseudo ucsc format
		island_libName1 = (peakfile).split('/')[-1]
		island_suffix1 = island_libName1.split('.')[-1] 
		island_libName1 = island_libName1.split('.')[0]
		island_extension1 = "-" + island_libName1 + '.' + island_suffix1 + "1"
		SeparateByChrom.separateByChrom(chroms, peakfile, island_extension1)
	else:
		print peakfile, " is not found";
		sys.exit(1)
	
	for chrom in chroms: 
		if chrom in entrez_genes.chroms:
			entrez_genes_by_chrom =  Entrez.KnownEntrezGenes([chrom], entrez_genes.subset_by_chrom(chrom))
			this_chrom_length = chrom_lengths[chrom]
			
			# Load in the PA peak information 
			if Utility_extended.fileExists(chrom + island_extension1):
				inf = open(chrom + island_extension1, 'r')
				# Read in the peaks and separate the forward strand peaks and the reverse strand peaks
				five_peaks = [] # peaks on forward strand, element (location, read_count)
				three_peaks = [] # peaks on reverse strand, element (location, read_count)
				for line in inf:
					line = line.strip();
					sline = line.split();
					strand = sline[2]
					if plus.match(strand):
						if float(sline[10]) >= peak_threshold:
							five_peaks.append ((int(sline[3]), float(sline[10])))
					elif minus.match(strand):
						if float(sline[10]) >= peak_threshold:
							three_peaks.append ((int(sline[4]), float(sline[10])))
				five_peaks = sorted(five_peaks, key = itemgetter(0)) #sort according to location
				five_peaks_location = [item[0] for item in five_peaks]
				three_peaks = sorted(three_peaks, key = itemgetter(0))
				three_peaks_location = [item[0] for item in three_peaks]
				inf.close()
				
				for entrez_id in entrez_genes_by_chrom.entrez_ids:
					gene = entrez_genes_by_chrom.entrez_genes[entrez_id] # an EntrezGene class object
					
					# For the set of transcripts, use the longest 3UTR at the designated representative 3UTR
					transcript_with_longest_3UTR = gene.identify_transcript_with_longest_3UTR() # a UCSC class object
					
					if plus.match(transcript_with_longest_3UTR.strand):
						start = transcript_with_longest_3UTR.cdsEnd
						end = min(transcript_with_longest_3UTR.txEnd + downstream_extension, this_chrom_length)
						start_ind = bisect.bisect_left(five_peaks_location, start);
						end_ind = bisect.bisect_right(five_peaks_location, end);
						peaks_on_3UTR = five_peaks[start_ind: end_ind] #[(mode_location, readcount)]
					if minus.match(transcript_with_longest_3UTR.strand):
						start = max(transcript_with_longest_3UTR.txStart - downstream_extension, 0)
						end = transcript_with_longest_3UTR.cdsStart
						start_ind = bisect.bisect_left(three_peaks_location, start);
						end_ind = bisect.bisect_right(three_peaks_location, end);
						peaks_on_3UTR = three_peaks[start_ind: end_ind]
					ThreeUTR_length = end - start + 1 #length includes the downstream extension
					peaks_on_entrez_3UTRs[entrez_id] = (gene, ThreeUTR_length, peaks_on_3UTR)
				
	SeparateByChrom.cleanup(chroms, island_extension1)
	return peaks_on_entrez_3UTRs
						
						
def Calculate3UTRUsage(peaks_on_entrez_3UTRs, outfile):
	
	"""
	Calculate the characteristics of 3UTR usage by mean, standard deviation and totoal read count, from the peaks  
	"""
	zero_pa_genes = 0
	
	outf = open(outfile, 'w')
	#outline = "# Refseq ID is used to denote the 3UTR used"
	#outf.write(outline)
	#outline = "# Entrez ID \t Main Refseq ID \t 3UTR length(inc. ext) \t Length Index \t PA Multiplicity Index \t Peaks Read Count \t RefSeq IDs \t Gene symbols \n"
	#outf.write(outline)
	
	for entrez_id in peaks_on_entrez_3UTRs.keys():
		gene = peaks_on_entrez_3UTRs[entrez_id][0]
		ThreeUTR_length = peaks_on_entrez_3UTRs[entrez_id][1]
		peaks = peaks_on_entrez_3UTRs[entrez_id][2]
		
		results = Characteristics_by_gene(gene, peaks) #(total_read_count, mean, pmi)
		if results[0] == 0:
			zero_pa_genes += 1
			
		gene_symbol = []
		for mygene in genes:
			if mygene.additional_annotations[0] not in gene_symbol:
				gene_symbol.append(mygene.additional_annotations[0])
		outline = str(entrez_id) + "\t" + gene.name + "\t" + str(ThreeUTR_length) + "\t" + "\t".join(map(str, results[1:])) + "\t" + str(results[0])  + "\t" + ",".join([gene.name for gene in genes]) + "\t" + ','.join(gene_symbol) + "\n"
		#print outline
		outf.write(outline)
			
			#if entrez_id == 182:
				#print gene.getAll()
				#print
				#for mygene in genes:
					#print mygene.getAll()
				#print
				#print three_peaks[start_ind: end_ind]
				#print 
				#print gene.txEnd
				#print results
	print "Number of entrez genes with no PA reads on longest 3UTR: ", zero_pa_genes				


def Characteristics_by_gene (gene, listofpeaks):
	# gene is a entrez_gene object
	# list of peaks (position, read_count)
	# can also be used for read-based calculation, where each peak is represented by (position, 1)
	# output: total read count, mean distance away from annotated PA, SD 
	transcript_with_longest_3UTR = gene.identify_transcript_with_longest_3UTR()
	if plus.match(gene.strand):
		annotated_PAS = transcript_with_longest_3UTR.txEnd
	elif minus.match(gene.strand):
		annotated_PAS = transcript_with_longest_3UTR.txStart
	return Characteristics(annotated_PAS, gene.strand, listofpeaks)
	
def Characteristics(annotated_PAS, strand, listofpeaks):
	"""
	annotated_PAS: a location
	# list of peaks (position, read_count)
	# can also be used for read-based calculation, where each peak is represented by (position, 1)
	# output: total read count, mean distance away from annotated PA, SD 
	"""
	
	if len(listofpeaks) == 0:
		return (0.0, 0.0, 0.0)
	else:
		read_counts = [item[1] for item in listofpeaks]
		total_read_count = float(sum(read_counts))
		if plus.match(strand):
			mean = sum ([(annotated_PAS - item[0] )*item[1]/total_read_count for item in listofpeaks])
			pmi = sum([(item[0] - (annotated_PAS - mean))*(item[0] - (annotated_PAS - mean))*item[1]/total_read_count for item in listofpeaks])
		elif minus.match(strand):
			mean = sum([( item[0] - annotated_PAS )*item[1]/total_read_count for item in listofpeaks])
			pmi = sum([(item[0] - (annotated_PAS + mean))*(item[0] - (annotated_PAS + mean))*item[1]/total_read_count for item in listofpeaks])
		pmi = sqrt(pmi)
		return (total_read_count, mean, pmi)

		
		
def Characteristics_by_Tian(gene, listofpeaks, allowance=50):
	"""
	
	The idea is to 
	
	# gene is a entrez_gene object
	# list of peaks (position, read_count)
	
	needs additional thinking? 
	
	Use the annotated PA, and the rest? Or the Biggest peak and the rest?
	How about those with elongations?
	
	The module implement the approach used  by Tian and Burge, which basically only counts two parts: 1) annotated PA, and 2) the rest
	"""
	
	if plus.match(gene.strand):
		annotated_PAS = gene.txEnd
	elif minus.match(gene.strand):
		annotated_PAS = gene.txStart


def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--peakfile", action="store", type="string", dest="peakfile", help="input ucsc file for PA peaks ", metavar="<file>")
	parser.add_option("-u", "--annotationfile", action="store", type="string", dest="annotationfile", help="pickle file for annotations ", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="outfile name", metavar="<file>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-t", "--peak_threshold", action="store", type="int", dest="peak_threshold",help="Peak threshold", metavar="<int>")
	parser.add_option("-d", "--3UTRdownstreamextension", action="store", type="int", dest="downstream_extension",help="3UTR down stream extension", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
		parser.print_help()
		sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	# entrez_gene_collection is a KnownEntrezGenes class object. The core is a entrez_genes.entrez_genes, a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
	
	# test module
	test = 0
	if test == 1:
		print "Testing gene structure"
		test_id = 54
		Entrez.test_gene_structure(entrez_gene_collection, test_id)

	# Filter cluster of refseq_ids (keyed by entrez_id) according to the criterion of identical cdsEnd
	entrez_ids_with_unique_cdsEnd = entrez_gene_collection.get_ids_unique_cdsEnd()
	print "There are ", len(entrez_ids_with_unique_cdsEnd), " Entrez IDs each of which has a unique cdsEnd."
	
	# Additional filter to remove clusters with intron-containing 3UTRs
	allowance=0
	ids=entrez_ids_with_unique_cdsEnd
	entrez_ids_with_intronless_3UTRs = entrez_gene_collection.get_ids_with_intronless_3UTR(allowance, ids)
	print "There are %d Entrez_ids with additional requirement of intronless 3UTR: ", %(len(entrez_ids_with_intronless_3UTRs))
	
	entrez_gene_subset = Entrez.KnownEntrezGenes(chroms, entrez_gene_collection.subset(entrez_ids_with_intronless_3UTRs))
	
	peaks_on_entrez_3UTRs = AssignPeaksToEntrez3UTRs(entrez_gene_subset, opt.peakfile, chroms, chrom_lengths, opt.peak_threshold, opt.downstream_extension)
	
	output = open(libName + "_PA_Peaks_associated_with_Annotations.pkl", 'wb')
	pickle.dump(peaks_on_entrez_3UTRs, output)
	output.close()
	
	Calculate3UTRUsage(peaks_on_entrez_3UTRs, final_entrez_id_collection, opt.outfile)
	
	
if __name__ == "__main__":
	main(sys.argv)