import HTSeq
import sys, os, re
from string import *
from optparse import OptionParser
import collections
import itertools

def get_read_length(alt_iv_seq):
	read_length = 0
	for alt_iv in alt_iv_seq:
		read_length += alt_iv.length
	return read_length

def reverse_strand( iv ):
	if iv.strand == "+":
		iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
	elif iv.strand == "-":
		iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
	return iv_reversed

def combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq):
	combine_alt_iv_seq = []
	alt_iv_seq = alt_first_iv_seq + alt_second_iv_seq
	alt_iv_seq.sort(key=lambda x: x.start)
	combine_alt_iv_seq.append(alt_iv_seq[0].copy())
	for alt in itertools.islice(alt_iv_seq, 1, None):
		if alt.overlaps(combine_alt_iv_seq[-1]):
			combine_alt_iv_seq[-1].extend_to_include(alt)
		else:
			combine_alt_iv_seq.append(alt.copy())
	return combine_alt_iv_seq 

def iv_to_str(iv):
	return iv.chrom + ':' + str(iv.start) + '-' + str(iv.end)   

def main(argv):
	parser = OptionParser()
	
	parser.add_option("-b", "--bamfile", action="store", type="string", dest="bamfile", metavar="<file>", help="RNA seq read file")
	parser.add_option("-p", "--read-type", action="store", type="string", dest="readtype", metavar="<file>", help="read type: single or paired", default='single')
	parser.add_option("-s", "--library-type", action="store", type="string", dest="libtype", metavar="<file>", help="Library type. DEFAULT: \"fr-unstranded\" (unstranded). Use \"fr-firststrand\" or \"fr-secondstrand\" for strand-specific data", default='fr-unstranded')
	parser.add_option("-g", "--genes_GTF_file", action="store", type="string", dest="genesfile", metavar="<file>", help="GTF genes annotation file")
	parser.add_option("-k", "--gene_list", action="store", type="string", dest="gene_list", metavar="<file>", help="gene list")
	parser.add_option("-d", "--DOG_extension", action="store", type="int", dest="DOG_extension", help="downstream of gene extension",  metavar="<int>")	
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
		parser.print_help()
		sys.exit(1)
		
	if opt.gene_list:	
		gene_list = []	
		f = open(opt.gene_list, "r")
		for line in f:
			line = line.strip()
			sline = line.split()
			gene_list.append(sline[0])
			
	if opt.libtype == "fr-firststrand" or opt.libtype == "fr-secondstrand":
		genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)
		DOGs = HTSeq.GenomicArrayOfSets("auto", stranded=True)
	elif opt.libtype == "fr-unstranded":
		genes = HTSeq.GenomicArrayOfSets("auto", stranded=False)
		DOGs = HTSeq.GenomicArrayOfSets("auto", stranded=False)
		
	DOG_extension = opt.DOG_extension
	
	counts = collections.Counter()
	DOG_length_dict = collections.Counter()
	DOG_iv_str_dict = collections.defaultdict(lambda : '')
	total_count = 0
	num_transcripts_dict = {}

	
	if opt.gene_list:
		gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
		for feature in gtffile:
			if feature.type == 'gene_region' and feature.attr["gene_id"] in gene_list:
				gene_region_iv = feature.iv
				num_transcripts_dict[feature.attr["gene_id"]] = int(feature.attr["transcripts_in_gene"])
				genes[gene_region_iv] += 'gene_region'
				
			elif feature.type == 'transcript_region' and feature.attr["gene_id"] in gene_list:
				transcript_region_iv = feature.iv
				if transcript_region_iv.strand == "+":
					DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.end, transcript_region_iv.end + DOG_extension, '+')
				elif transcript_region_iv.strand == "-":
					if transcript_region_iv.start - DOG_extension < 0:
						DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.start, '-')
					else:
						DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start - DOG_extension, transcript_region_iv.start, '-')

				gene_id = feature.attr["gene_id"]
				transcript_id = feature.attr["transcript_id"]
				
				genes[DOG_iv] += gene_id	
					
	else:
		gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
		for feature in gtffile:
			if feature.type == 'gene_region':
				gene_region_iv = feature.iv
				num_transcripts_dict[feature.attr["gene_id"]] = int(feature.attr["transcripts_in_gene"])
				genes[gene_region_iv] += 'gene_region'
				
			elif feature.type == 'transcript_region':
				transcript_region_iv = feature.iv
				if transcript_region_iv.strand == "+":
					DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.end, transcript_region_iv.end + DOG_extension, '+')
				elif transcript_region_iv.strand == "-":
					if transcript_region_iv.start - DOG_extension < 0:
						DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.start, '-')
					else:
						DOG_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start - DOG_extension, transcript_region_iv.start, '-')

				gene_id = feature.attr["gene_id"]
				transcript_id = feature.attr["transcript_id"]
				
				genes[DOG_iv] += gene_id		
		
			
	for iv, step_set in genes.steps():
		if len(step_set) != 0 and 'gene_region' not in step_set:
			for gene_id in list(step_set):
				DOGs[iv] += gene_id
				DOG_length_dict[gene_id] += iv.length
				DOG_iv_str_dict[gene_id] += iv_to_str(iv) + ','	
				
				if gene_id not in counts.keys():
					counts[gene_id] = 0									
				
	bamfile = HTSeq.BAM_Reader(opt.bamfile)
	
	# Single end
	if opt.readtype == "single":	
		for alt in bamfile:
			if alt.aligned and alt.optional_field("NH") == 1 and re.match('chr', alt.iv.chrom):
				total_count += 1
				
				if opt.libtype == "fr-secondstrand" or opt.libtype == "fr-unstranded":
					alt_iv_seq = [ co.ref_iv for co in alt.cigar if co.type == "M" and co.size > 0 ]
				elif opt.libtype == "fr-firststrand":
					alt_iv_seq = [ reverse_strand(co.ref_iv) for co in alt.cigar if co.type == "M" and co.size > 0 ]
				
				read_length = get_read_length(alt_iv_seq)
				for alt_iv in alt_iv_seq:
					for iv, step_set in DOGs[alt_iv].steps():
						if step_set:
							num_genes = len(step_set)
							for gene_id in list(step_set):
								counts[gene_id] += iv.length * 1.0 / read_length / num_genes
								
	elif opt.readtype == "paired":
		for alt_first, alt_second in HTSeq.pair_SAM_alignments(bamfile):
			if alt_first == None or alt_second == None:
				continue
			if alt_first.aligned and alt_first.optional_field("NH") == 1 and alt_second.aligned and alt_second.optional_field("NH") == 1 and re.match('chr', alt_first.iv.chrom) and re.match('chr', alt_second.iv.chrom):
				if alt_first.iv.chrom != alt_second.iv.chrom:
					continue

				total_count += 1

				# If strand_specific is "yes", first mate maps to the transcript strand, and the second mate maps to the opposite strand
				# Reverse the second mate to the transcript strand                                        
				if opt.libtype == "fr-secondstrand" or opt.libtype == "fr-unstranded":
					alt_first_iv_seq = [ co.ref_iv for co in alt_first.cigar if co.type == "M" and co.size > 0 ]
					alt_second_iv_seq = [ reverse_strand(co.ref_iv) for co in alt_second.cigar if co.type == "M" and co.size > 0 ]

				# If strand_specific is "reverse", first mate maps to the opposite strand, and the second mate maps to the transcript strand
				# Reverse the first mate to the transcript strand                                                 
				elif opt.libtype == "fr-firststrand":
					alt_first_iv_seq = [ reverse_strand(co.ref_iv) for co in alt_first.cigar if co.type == "M" and co.size > 0 ]
					alt_second_iv_seq = [ co.ref_iv for co in alt_second.cigar if co.type == "M" and co.size > 0 ]
				alt_iv_seq = combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)	
				
				read_length = get_read_length(alt_iv_seq)
				for alt_iv in alt_iv_seq:
					for iv, step_set in DOGs[alt_iv].steps():
						if step_set:
							num_genes = len(step_set)
							for gene_id in list(step_set):
								counts[gene_id] += iv.length * 1.0 / read_length / num_genes				
				
				
	f = open(opt.outfile, "w")
	outline = "gene_id" +  '\t' + "DOG_iv" + "\t" + "DOG_length" + "\t" + "DOG_read_count" + '\t' + "RPKM" + '\n'
	f.write(outline)
	
	for gene_id in sorted(counts.keys()):
		RPKM = counts[gene_id] / (DOG_length_dict[gene_id] / 1000.0) / (total_count / 1000000.0)
		outline = gene_id + "\t" + DOG_iv_str_dict[gene_id] + "\t" + str(DOG_length_dict[gene_id]) + "\t"+ str(counts[gene_id]) + "\t" + str(RPKM) + "\n"
		f.write(outline)
		
	f.close()			
		
if __name__ == "__main__":
	main(sys.argv)
					
		
		
				
		