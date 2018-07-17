#!/usr/bin/env python

# 
# class version of the Generate_Entrez_UCSC.py 
# cluster the refseq ids to Entrez id. 
# Remove collision (to reduce RNA contamination in intronic regions) 
# Remove NR genes (if any refseq_ID in the cluster is a NR)  and Non concordant clusters
# can also be used for clustering the refseq ids to gene symbol.

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import time
import bisect
from operator import itemgetter
import copy
try:
   import cPickle as pickle
except:
   import pickle

sys.path.append('/home/data/SICER1.1/SICER/lib')
sys.path.append('/home/data/SICER1.1/SICER/extra')

import GenomeData 
import UCSC_revised
import Utility_extended

plus = re.compile("\+");
minus = re.compile("\-");
comment = re.compile("#|track")
NR = re.compile("NR|nr")

extension = 5000

class EntrezGene:
	"""
	self.transcripts: [transcript], transcript is an instance of a UCSC class object
	exonic region: internal exons.
	"""
	#def __init__(self, ucsc_transcripts, extension = 0):
	def __init__(self, ucsc_transcripts):
		self.transcripts = copy.deepcopy(ucsc_transcripts)
		self.num_transcripts = len(ucsc_transcripts)
		
		self.concordance = self.check_concordance(extension) 
		if self.concordance == 1:
			self.chrom = self.get_chrom()
			self.strand = self.get_strand()
			self.boundaries = UCSC_revised.getBoundariesMergedTranscripts(self.transcripts)# [(start, end)], start < end
			self.merged_exonic_regions = UCSC_revised.getMergedExonicRegions(self.transcripts)# Return the merged exons in the format of list of [(start, end)]
			self.shared_exonic_regions = UCSC_revised.getSharedExonicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_intronic_regions = UCSC_revised.getSharedIntronicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_3UTR = UCSC_revised.getShared3UTR(self.transcripts, 5) # minimum_width=5
			self.shared_5UTR = UCSC_revised.getShared5UTR(self.transcripts, 5) # minimum_width=5
	
	def insert(self, ucsc_gene):
		self.transcripts.append(ucsc_gene)
		self.num_transcripts = len(self.transcripts)
		self.concordance = self.check_concordance(extension) 
		if self.concordance == 1:
			self.chrom = self.get_chrom()
			self.strand = self.get_strand()
			self.boundaries = UCSC_revised.getBoundariesMergedTranscripts(self.transcripts)# [(start, end)]
			self.merged_exonic_regions = UCSC_revised.getMergedExonicRegions(self.transcripts)# Return the merged exons in the format of list of [(start, end)]
			self.shared_exonic_regions = UCSC_revised.getSharedExonicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_intronic_regions = UCSC_revised.getSharedIntronicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_3UTR = UCSC_revised.getShared3UTR(self.transcripts, 5) # minimum_width=5
			self.shared_5UTR = UCSC_revised.getShared5UTR(self.transcripts, 5) # minimum_width=5
			
	def remove(self, ucsc_gene):
		self.transcripts.remove(ucsc_gene)
		self.num_transcripts = len(self.transcripts)
		self.concordance = self.check_concordance(extension) 
		if self.concordance == 1:
			self.chrom = self.get_chrom()
			self.strand = self.get_strand()
			self.boundaries = UCSC_revised.getBoundariesMergedTranscripts(self.transcripts)# [(start, end)]
			self.merged_exonic_regions = UCSC_revised.getMergedExonicRegions(self.transcripts)# Return the merged exons in the format of list of [(start, end)]
			self.shared_exonic_regions = UCSC_revised.getSharedExonicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_intronic_regions = UCSC_revised.getSharedIntronicRegions(self.transcripts, 5) # minimum_width=5
			self.shared_3UTR = UCSC_revised.getShared3UTR(self.transcripts, 5) # minimum_width=5
			self.shared_5UTR = UCSC_revised.getShared5UTR(self.transcripts, 5) # minimum_width=5
		
	def check_concordance(self, my_extension=0):
		"""
		Check the concordance of this list of transcripts 
		extension is on both sides
		"""
		if self.check_chrom_concordance() == 0:
			return 0
		elif self.check_strand_concordance() == 0: #already considers requiring the chrom concordance
			return 0
		elif self.check_coordinate_concordance(my_extension) == 0:
			return 0
		else:
			return 1
	
	def eligibility(self, my_extension=0):
		"""
		Determines the eligilibity of a gene
		
		Currently:
		concordant: these transcripts need to have overlap of some sorts, have to be on the same chrom, same strand, 
		NR: non-coding genes are excluded
		
		"""
		if self.check_concordance(my_extension) == 1:
			if self.is_NR_containing_gene() == 0:
				return 1
			else: 
				return 0
		else:
			return 0
	
	def check_strand_concordance(self):
		"""
		If transcripts with the same Entrez ID are not on the same strand, they are not concordant
		"""
		mystrand = self.transcripts[0].strand
		if self.num_transcripts > 1:
			for transcript in self.transcripts:
				if transcript.strand != mystrand:
					return 0
		return 1

	def check_coordinate_concordance(self, my_extension):
		"""
		If any two extended transcripts with the same Entrez ID have zero overlaps, they are not concordant. 
		"""
		if self.num_transcripts > 1:
			for index_1 in xrange(self.num_transcripts-1):
				for index_2 in range(index_1+1, self.num_transcripts):
					if overlap(self.transcripts[index_1], self.transcripts[index_2], my_extension) == 0:
						return 0
		return 1

	def check_chrom_concordance(self):
		chrom = self.transcripts[0].chrom
		if self.num_transcripts > 1:
			for transcript in self.transcripts:
				if transcript.chrom != chrom:
					return 0
		return 1
	
	def is_NR_containing_gene(self):
		"""
		This is non-trivial, there are cases that under the same entrez_ID, one is NR and the other is NM To be conservative, let us discard the cluster if any of the refseq ID is NR.
		"""
		mynr = 0
		for transcript in self.transcripts:
			if NR.match(transcript.name):
				mynr = 1
		return mynr
	
	def is_NR_gene(self):
		"""
		Find Genes that includes only NR transcripts
		"""
		mynr = 1 
		for transcript in self.transcripts:
			if not NR.match(transcript.name):
				mynr = 0
		return mynr
		
	def get_chrom(self):
		if self.check_chrom_concordance() == 1:
			chrom = self.transcripts[0].chrom
			return chrom 
		else:
			print "Transcripts are on different chroms"
			return ""
			
	def get_number_of_exons(self):
		"""
		return the number of internal exons: [number of exons]
		"""
		number_of_exons = []
		for mytranscript in self.transcripts:
			number_of_exons.append(len(mytranscript.getExons()))
		return number_of_exons
	
	def get_transcript_lengths(self):
		"""
		return the transcript lengths: [transcript length]
		"""
		transcript_lengths=[]
		upstream = 0
		downstream = 0
		for mytranscript in self.transcripts:
			(egb_start, egb_end) = mytranscript.getExtendedGeneBody(upstream, downstream)
			transcript_lengths.append( abs(egb_start -egb_end) )
		return transcript_lengths
	
	def get_strand(self):
		return self.transcripts[0].strand
		
	def info(self):
		print "There are %d transcripts" % self.num_transcripts
		if self.concordance == 1:
			print "Chrom: ", self.chrom 
			print "Strand: ", self.strand 
			print "Boundaries: ", self.boundaries 
			print "Merged exonic regions: ", self.merged_exonic_regions 
			print "Shared_exonic_regions: ", self.shared_exonic_regions # minimum_width=5
			print "Shared_intronic_regions: ", self.shared_intronic_regions  # minimum_width=5
			print "Shared 3UTR: ", self.shared_3UTR  # minimum_width=5
			print "Shared 5UTR: ", self.shared_5UTR  # minimum_width=5
		for transcript in self.transcripts:
			print transcript.getAll()
	
	def isIdenticalCdsEnd(self, allowance=10):
		"""
		Determines whether this cluster of transcript share the same cdsEnd
		"""
		answer = 1 
		if plus.match(self.strand):
			cdsends = [ item.cdsEnd for item in self.transcripts]	
		elif minus.match(self.strand):
			cdsends = [ item.cdsStart for item in self.transcripts]
		mean = sum(cdsends) / (len(cdsends) * 1.0)
		max_dist = max([abs(item-mean) for item in cdsends])
		if max_dist > allowance:
			answer = 0
		return answer
	
	def get_3UTRs(self, downstream = 0):
		"""
		Returns 3UTRs as a list of tuples. 
		For a gene in positive strand, [(cdsend, txend)]
		For a gene in negative strand  [(txstart, cdsStart)]
		
		"""
		three_UTRs = []
		for transcript in self.transcripts:
			upstream = 0
			three_UTR = transcript.get3UTR(upstream, downstream)
			if  len(three_UTR) > 0:
				three_UTRs.append( three_UTR )
		return three_UTRs
	
	def identify_transcript_with_longest_3UTR(self):
		"""
		pick the longest 3UTR 
		return the corresponding transcript
		
		"""
		mymax = 0
		myindex = 0
		for index in xrange(len(self.transcripts)):
			transcript = self.transcripts[index]
			upstream = 0
			downstream = 0
			three_UTR = transcript.get3UTR(upstream, downstream) #(tUTR_start, tUTR_end)
			if len(three_UTR) > 0:
				tUTR_end = three_UTR[1]
				tUTR_start = three_UTR[0]
				length = tUTR_end - tUTR_start + 1
				if length > mymax:
					mymax = length
					myindex = index
		return copy.deepcopy(self.transcripts[myindex])
		
	def is_3UTR_intronic(self, allowance):
		"""
		allowance is for what?
		"""
		for transcript in self.transcripts:
			if transcript.Is3UTRIntronContaining(allowance) == 1:
				return 1
		return 0
		
		
class KnownEntrezGenes:
	"""Entrez genes: the core is {entrez_id:EntrezGene}"""
	def __init__(self, chroms, entrez_genes):
		self.entrez_genes = copy.deepcopy(entrez_genes) #dic keyed by entrez_id
		self.chroms = chroms
		self.entrez_ids = self.entrez_genes.keys() # list of entrez_ids
		
		self.num_genes = len((entrez_genes).keys())
		self.num_transcripts = 0
		for myid in self.entrez_ids:
			self.num_transcripts += (entrez_genes[myid]).num_transcripts
		print "There are %d entrez genes for %d refseq transcripts." % (self.num_genes, self.num_transcripts) 

		# {chrom:[entrez_ids]}
		self.entrez_ids_by_chrom = self.separate_ids_by_chrom()
		#{entrez_id:[id which is colliding with entrez_id]}
		extension = 10
		self.entrez_id_collision_dic = self.find_all_colliding_ids(extension)
		# gene_boundaries_with_ranking: {chrom:{(entrez_ids, position):ranking}}
		# ranked_gene_boundaries: {chrom:[(entrez_ids, position)]}
		self.gene_boundaries_with_ranking, self.ranked_gene_boundaries  = self.order_gene_boundaries()
		
	@classmethod
	def initiate_from_file(cls, chroms, file):
		"""
		Requires that the entrez_id to be the second of the gene.additional_annotations[1]!!
		"""
		refseq_genes = UCSC_revised.KnownGenes(file);
		print 'There are %d transcripts in %s' % (refseq_genes.getNumGenes(), file)
	
		total = 0
		total_total = 0
		non_concordant = 0
		entrez_collection = {}
		entrez_genes = {}
		
		#refseq ids in weird chroms are ignored
		mychroms = list( set(chroms) & set(refseq_genes.keys()) )
		
		# cluster refseq ids according to their Entrez id. 
		for chrom in mychroms:
			for gene in refseq_genes[chrom]:
				entrez_id = int(gene.additional_annotations[1]) #Note the column index!!!!!!!!!!!!!!
				if entrez_id not in entrez_collection.keys():
					entrez_collection[entrez_id]=[]
				entrez_collection[entrez_id].append(gene)
		total = len((entrez_collection).keys())
		for id in (entrez_collection).keys():
			total_total += len((entrez_collection)[id]) 
		print "After removing transcripts on exotic chroms, there are ", total, " Entrez IDs, covering ", total_total," transcripts. "
		
		for myid in entrez_collection.keys():
			entrez_gene = EntrezGene(entrez_collection[myid])
			#check concordance, which requires the transcripts to be on the same strand, and all pairs overlapping after extension
			if entrez_gene.concordance == 1: 
				entrez_genes[myid] = entrez_gene
			else:
				non_concordant += 1
		print "After removing non-concordant transcripts, there are ", total - non_concordant, " Entrez genes. "
		return cls(chroms, entrez_genes)


	def get_entrez_ids(self):
		"""
		returns the list of entrez ids in entrez_genes
		"""
		return self.entrez_ids
	
	def get_total_entrez_ids(self):
		return self.num_genes
		
	def get_num_transcript(self):
		return self.num_transcripts
	
	def check_id_set(self, ids=None):
		"""
		check whether [ids] is a subset of self.entrez_genes
		if not, return intersection
		"""
		if ids is None:
			universe = self.entrez_ids
		else:
			universe = list(set(ids) & set(self.entrez_ids))
			if len(universe) < len(ids):
				print "%d ids are thrown out" % ( len(ids) - len(universe))
		return universe
	
	def filter_by_eligibility(self, ids=None, extension=0):
		"""
		Output the list of eligible ids
		"""
		myids = self.check_id_set(ids)
		
		eligible_ids = []
		for myid in myids:
			entrez_gene = (self.entrez_genes[myid])
			if entrez_gene.eligibility(extension) == 1:
				eligible_ids.append(myid)
		print "The number of ids after filtering by eligibility (including corcordance check and NR removal)is ", len(eligible_ids)
		return eligible_ids
		
	def output(self, outfile, ids=None):
		"""
		Output entrez genes represented by [ids] in file
		"""
		myids = self.check_id_set(ids)
		myids = sorted(myids)
		
		outf = open(outfile, 'w')
		for entrez_id in myids:
			transcripts = (self.entrez_genes[entrez_id]).transcripts #[ucsc]
			for transcript in transcripts:
				outline = transcript.getAll()
				outf.write(outline)
		print "The number of genes output into %s is %d"  % (outfile, len(myids))
		outf.close()
		
	def output_pickle(self, outfilename, ids=None):
		"""
		pickle entrez genes represented by {id:entrez_gene} in file
		"""
		myids = self.check_id_set(ids)
		mysubset = self.subset(myids)
		output = open(outfilename, 'wb')
		pickle.dump(mysubset, output)
		print "The number of genes output to %s is %d " % (outfilename, len(myids))
		output.close()	
		
	def subset(self, ids):
		# ids is a list of entrez_ids
		# return a entrez_genes dic
		myset=self.extract(ids, self.entrez_genes)
		return myset
	
	def get_transcript_ids(self, entrez_ids=None):
		myids = self.check_id_set(entrez_ids)
		transcript_id_set = []
		for myid in myids:
			gene = self.entrez_genes[myid]
			names = [transcript.name for transcript in gene.transcripts]
			transcript_id_set.extend(names)
		return transcript_id_set
		
	def separate_ids_by_chrom(self, ids=None):
		"""
		Generate a dic {chrom:[entrez_ids]}
		"""
		myids = self.check_id_set(ids)
		entrez_ids_by_chrom = {}
		for entrez_id in myids:
			gene = self.entrez_genes[entrez_id]
			if gene.chrom in self.chroms:
				if gene.chrom not in entrez_ids_by_chrom.keys():
					entrez_ids_by_chrom[gene.chrom] = []
				entrez_ids_by_chrom[gene.chrom].append(entrez_id)
		return entrez_ids_by_chrom
	
	def subset_by_chrom(self, chrom, ids=None):
		# myids is a list of entrez_ids 
		# return a entrez_genes dic: {entrez_id:entrez_gene}
		if chrom in self.chroms:
			myids = (self.separate_ids_by_chrom(ids))[chrom]	
			return self.subset(myids)
		else:
			return {}

	def id_clustering_histogram(self, outfile):
		"""
		histogram
		"""
		outf = open(outfile, 'w')
		setmax = 0
		for myid in (self.entrez_genes).keys():
			setmax = max (setmax, (self.entrez_genes[myid]).num_transcripts)
		entrez_id_collection_histogram = [0]*(setmax+1)
		for myid in self.entrez_genes.keys():
				entrez_id_collection_histogram[(self.entrez_genes[myid]).num_transcripts] += 1
		for index in xrange(len(entrez_id_collection_histogram)):
			outline = str(index) + '\t' + str(entrez_id_collection_histogram[index]) + '\n'
			outf.write(outline)
		outf.close()

	def find_simple_genes(self, ids=None):
		"""
		Genes with only one annotated transcripts
		return [ids]
		"""
		myids = self.check_id_set(ids)
		simple_gene_ids= []
		for entrez_id in myids:
			entrez = self.entrez_genes[entrez_id]
			if len(entrez.transcripts) == 1:
				simple_gene_ids.append(entrez_id)
		return simple_gene_ids		

	def get_strand_specific_ids(self, strand, ids=None):
		"""
		returns [ids]
		"""
		myids = self.check_id_set(ids)
		mygenes=[]
		assert (strand == "+" or strand =="-")
		for entrez_id in myids:
			entrez = self.entrez_genes[entrez_id]
			if entrez.strand == strand:
				mygenes.append(entrez_id)
		return mygenes

	def extract(self, keys, d):
		return dict((k, d[k]) for k in keys if k in d)
	
	def get_ids_with_unique_cdsEnd(self, allowance=10, ids=None):
		"""
		Extract those entrez genes whose transcripts share the same cds END
		For PA analysis
		returns [ids]
		"""
		myids = self.check_id_set(ids)
		final_ids = []
		for entrez_id in myids:
			gene = self.entrez_genes[entrez_id]
			if gene.isIdenticalCdsEnd(allowance) == 1:
				final_ids.append(entrez_id)
		return final_ids
	
	def get_ids_with_intronless_3UTR(self, allowance=0, ids=None):
		"""
		Extract only those genes with intronless_3UTR. 
		"""
		myids = self.check_id_set(ids)
		final_ids = []
		for entrez_id in myids:
			gene = self.entrez_genes[entrez_id]
			if gene.is_3UTR_intronic(allowance) == 0:
				final_ids.append(entrez_id)
		return final_ids
		
	def order_gene_boundaries_on_a_chrom (self, chrom, ids=None):
		
		"""
		arrange all boundaries positions in a sorted list of tuples, and store the association of tuple and ranking in a dictionary so that the ranking can be found by (entrez_id, location)
		
		returns 
		gene_boundaries_with_ranking: {(entrez_id, position): ranking along the genome}
		ranked_gene_boundaries: [(entrez_id, position)]
		"""
		myids = list(set(self.check_id_set(ids)) & set(self.entrez_ids_by_chrom[chrom]))
		tags =[] # [(entrez_id, location)]
		gene_boundaries_with_ranking = {}
		
		for entrez_id in myids:
			boundaries = (self.entrez_genes)[entrez_id].boundaries[0]
			start = boundaries[0]
			tags.append((entrez_id, start))
			end = boundaries[1]
			tags.append((entrez_id, end))
		ranked_gene_boundaries = sorted(tags, key=itemgetter(1))
		for index in xrange(len(ranked_gene_boundaries)):
			entrez_id = ranked_gene_boundaries[index][0]
			position = ranked_gene_boundaries[index][1] # Can be start or end
			gene_boundaries_with_ranking[(entrez_id, position)] = index
		return gene_boundaries_with_ranking, ranked_gene_boundaries
		
	def order_gene_boundaries(self, ids=None):
		"""
		return 
		gene_boundaries_with_ranking: {chrom:{(entrez_id, position): ranking along the genome}}
		ranked_gene_boundaries: {chrom:[(entrez_id, position)]}
		"""
		gene_boundaries_with_ranking = {}
		ranked_gene_boundaries = {}
		for chrom in self.chroms:
			gene_boundaries_with_ranking[chrom], ranked_gene_boundaries[chrom] = self.order_gene_boundaries_on_a_chrom(chrom, ids)
		return gene_boundaries_with_ranking, ranked_gene_boundaries
	
	def get_neighbor_ids_of_a_gene(self, entrez_id, number, max_distance = -1, extension = 10):
		"""
		caveat: the neighbors found are limited within the confines of chosen entrez_ids, those not picked in defining the universe of genes are neglected
		
		Those genes overlapping are not counted as neighbors 
		
		max_distance = -1 means no distance requirement
		
		return:
		{"upstream": [entrez_id]; "downstream": [entrez_id]}
		"""
		mygene = (self.entrez_genes)[entrez_id]
		mystart = mygene.boundaries[0][0] #mygene.boundaries: [(start, end)]
		myend = mygene.boundaries[0][1] # end always behind start
		
		chrom = mygene.chrom
		
		my_gene_boundaries_with_ranking = self.gene_boundaries_with_ranking[chrom] #{(entrez_id, boundary_position): ranking along the genome)}
		my_ranked_gene_boundaries = self.ranked_gene_boundaries[chrom] #[(entrez_id, boundary_position)]
		
		if plus.match(mygene.strand):
			
			fiveprime_ranking = my_gene_boundaries_with_ranking[(entrez_id, start)]
			neighbor = max(0, fiveprime_ranking - 1)
			upstream_ids = set([])
			while len(upstream_ids) < number:
				neighbor_entrez_id = my_ranked_gene_boundaries[neighbor][0]
				# Overlapping genes are removed
				if is_overlapping(neighbor_entrez_id, entrez_id, extension) == 0: # no overlap
					upstream_ids.add(neighbor_entrez_id)
				if neighbor == 0:
					break
				neighbor -= 1
					
			threeprime_ranking = my_gene_boundaries_with_ranking[(entrez_id, end)]
			neighbor = min(len(my_ranked_gene_boundaries) - 1, threeprime_ranking + 1)
			downstream_ids = set([])
			while len(downstream_ids) < number:
				neighbor_entrez_id = my_ranked_gene_boundaries[neighbor][0]
				# Overlapping genes are removed
				if is_overlapping(neighbor_entrez_id, entrez_id, extension) == 0: # no overlap
					downstream_ids.add(neighbor_entrez_id)
				if neighbor == len(my_ranked_gene_boundaries) - 1:
					break
				neighbor += 1
			
		elif minus.match(mygene.strand):
			fiveprime_ranking = my_gene_boundaries_with_ranking[(entrez_id, end)]
			neighbor = min(len(my_ranked_gene_boundaries) - 1, fiveprime_ranking + 1)
			upstream_ids = set([])
			while len(upstream_ids) < number:
				neighbor_entrez_id = my_ranked_gene_boundaries[neighbor][0]
				# Overlapping genes are removed
				if is_overlapping(neighbor_entrez_id, entrez_id, extension) == 0: # no overlap
					upstream_ids.add(neighbor_entrez_id)
				if neighbor == len(my_ranked_gene_boundaries) - 1:
					break
				neighbor += 1
			
			threeprime_ranking = my_gene_boundaries_with_ranking[(entrez_id, start)]
			neighbor = max(0, threeprime_ranking - 1)
			downstream_ids = set([])
			while len(downstream_ids) < number:
				neighbor_entrez_id = my_ranked_gene_boundaries[neighbor][0]
				# Overlapping genes are removed
				if is_overlapping(neighbor_entrez_id, entrez_id, extension) == 0: # no overlap
					downstream_ids.add(neighbor_entrez_id)
				if neighbor == 0:
					break
				neighbor -= 1
		else:
			print "Strand information wrong"
			exit(1)
		
		upstream_neighbour_ids = list(upstream_ids)
		filtered_upstream_neighbour_ids = []
		downstream_neighbour_ids = list(downstream_ids)
		filtered_downstream_neighbour_ids = []
		
		if max_distance > 0:
			for entrez_id in upstream_neighbour_ids:
				gene = (self.entrez_genes)[entrez_id]
				start = gene.boundaries[0][0]
				end = gene.boundaries[0][1]
				distance = min(abs(start-mystart), abs(end-mystart), abs(start-myend), abs(end-myend))
				if distance <= max_distance:
					filtered_upstream_neighbour_ids.append(entrez_id)
			upstream_neighbour_ids = [myid for myid in filtered_upstream_neighbour_ids]
			for entrez_id in downstream_neighbour_ids:
				gene = (self.entrez_genes)[entrez_id]
				start = gene.boundaries[0][0]
				end = gene.boundaries[0][1]
				distance = min(abs(start-mystart), abs(end-mystart), abs(start-myend), abs(end-myend))
				if distance <= max_distance:
					filtered_downstream_neighbour_ids.append(entrez_id)
			downstream_neighbour_ids = [myid for myid in filtered_downstream_neighbour_ids]
		result = {}
		result["upstream"] = upstream_neighbour_ids
		result["downstream"] = downstream_neighbour_ids
		return result

			
			
			
	def get_intergenic_region_of_a_gene(self, entrez_id, extension=0):
		"""
		return {"upstream":region; "downstream":region} 
		"""
		mygene = (self.entrez_genes)[entrez_id]
		mystart = mygene.boundaries[0][0]
		myend = mygene.boundaries[0][1] # end always behind start
		
		neighbor = get_neighbor_ids_of_a_gene (entrez_id, 1, max_distance = -1, extension = 10)
		upstream_id = neighbor["upstream"][0]
		downstream_id = neighbor["downstream"][0]
		
		colliding_id_list = self.entrez_id_collision_dic[entrez_id]
		if len(colliding_id_list) > 0: #There is colliding id.
			mystart_in = 0
			myend_in = 0
			for myid in colliding_id_list:
				colliding_gene = self.entrez_genes[myid]
				start = colliding_gene.boundaries[0][0]
				end = colliding_gene.boundaries[0][1]
				if mystart >= start and mystart <= end:
					mystart_in = 1
				if myend >= start and myend <= end:
					myend_in = 1
					
			if mystart_in == 1:
				if plus.match(mygene.strand):
					upstream_region = (mystart - 1, mystart - 1)
				elif minus.match(mygene.strand):
					downstream_region = (mystart - 1, mystart - 1)
			else:
				if plus.match(mygene.strand):
					upstream_gene = (self.entrez_genes)[upstream_id]
					region_start = max (upstream_gene.boundaries[0][0], upstream_gene.boundaries[0][1])
					upstream_region = (region_start + 1, mystart - 1)
				elif minus.match(mygene.strand):
					downstream_gene = (self.entrez_genes)[downstream_id]
					region_start = max (downstream_gene.boundaries[0][0], downstream_gene.boundaries[0][1])
					downstream_region = (region_start + 1, mystart - 1)
			
			if myend_in == 1:
				if plus.match(mygene.strand):
					downstream_region = (myend + 1, myend + 1)
				elif minus.match(mygene.strand):
					upstream_region = (myend + 1, myend + 1)
			else:
				if plus.match(mygene.strand):
					downstream_gene = (self.entrez_genes)[downstream_id]
					region_end = min (downstream_gene.boundaries[0][0], downstream_gene.boundaries[0][1])
					downstream_region = (myend + 1, region_end - 1)
					
				elif minus.match(mygene.strand):
					upstream_gene = (self.entrez_genes)[upstream_id]
					region_end = min (upstream_gene.boundaries[0][0], upstream_gene.boundaries[0][1])
					upstream_region = (myend + 1, region_end - 1)
			
		else:#There is no colliding id
			upstream_gene = (self.entrez_genes)[upstream_id]
			if plus.match(mygene.strand):
				region_start = max (upstream_gene.boundaries[0][0], upstream_gene.boundaries[0][1])
				upstream_region = (region_start + 1, mystart - 1)
			elif minus.match(mygene.strand):
				region_end = min (upstream_gene.boundaries[0][0], upstream_gene.boundaries[0][1])
				upstream_region = (myend + 1, region_end - 1)
			
			downstream_gene = (self.entrez_genes)[downstream_id]
			if plus.match(mygene.strand):
				region_end = min (downstream_gene.boundaries[0][0], downstream_gene.boundaries[0][1])
				downstream_region = (myend + 1, region_end - 1)
			elif minus.match(mygene.strand):
				region_start = max (downstream_gene.boundaries[0][0], downstream_gene.boundaries[0][1])
				downstream_region = (region_start + 1, mystart - 1)
		result = {}
		result["upstream"] = upstream_region
		result["downstream"] = downstream_region
		return result
		
	def is_overlapping(self, entrez_id_a, entrez_id_b, extension=0):
		"""
		return 1 or 0
		"""
		gene_a = (self.entrez_genes)[entrez_id_a]
		a_start = gene_a.boundaries[0][0] - extension
		a_end = gene_a.boundaries[0][1] + extension
		assert a_start <= a_end
		
		gene_b = (self.entrez_genes)[entrez_id_b]
		b_start = gene_b.boundaries[0][0] - extension
		b_end = gene_b.boundaries[0][1] + extension
		assert b_start <= b_end
		
		if  a_end < b_start or b_end < a_start: #Non-overlap
			return 0
		else: 
			return 1
	
	
	def find_colliding_ids(self, entrez_gene_boundaries, extension=0):
		"""
		entrez_gene_boundaries: [(start, end, entrez_id)]
		build and return a dic: {id:[ids_in_collision]}, where the value is the list of ids which collides with  
		""" 
		
		#union the boundaries to find clusters
		clusters_of_ids = Utility_extended.union_with_trace(entrez_gene_boundaries, extension) #Output {(start, end):[ entrez_gene_boundaries elements that contribute to that region]}
		
		entrez_id_collision_dic = {}
		for item in entrez_gene_boundaries:
			myid = item[2]
			entrez_id_collision_dic[myid] = []
			
		for region in clusters_of_ids.keys():
			ids = [item[2] for item in clusters_of_ids[region]]
			# although in the same union, a pair of ids need not to be directly overlapping
			# Make sure a pair is traversed only once
			for i in xrange(len(ids)):
				myid = ids[i]
				for j in range(i+1, len(ids)):
					the_other_id = ids[j]
					if self.is_overlapping(myid, the_other_id, extension) == 1:
						entrez_id_collision_dic[myid].append(the_other_id)
						entrez_id_collision_dic[the_other_id].append(myid)
		return 	entrez_id_collision_dic			
	
	def find_all_colliding_ids(self, extension):
		"""
		every id has an entry. If no colliding ids: id:[]
		return {id:[ids_in_collision]}
		"""
		entrez_id_collision_dic = {}
		for chrom in self.chroms:
			#convert to [(start, end, annotation)]:
			entrez_gene_boundaries = []
			for entrez_id in self.entrez_ids_by_chrom[chrom]:
				myboundaries = (self.entrez_genes)[entrez_id].boundaries 
				start = myboundaries[0][0] # the genomic start irrespective of the direction of the gene
				end = myboundaries[0][1] # the genomic end irrespective of the direction of the gene
				entrez_gene_boundaries.append((start, end, entrez_id))
			entrez_id_collision_dic.update(self.find_colliding_ids(entrez_gene_boundaries, extension))
		return entrez_id_collision_dic
	
def overlap(transcript1, transcript2, my_extension=0):
	"""
	transcript1 and transcript2 are ucsc class object
	"""
	if transcript1.strand != transcript2.strand:
		return 0
	else: # the strand does not matter
		start1 = transcript1.txStart - my_extension
		start2 = transcript2.txStart - my_extension
		end1 = transcript1.txEnd + my_extension
		end2 = transcript2.txEnd + my_extension
		assert start1 <= end1
		assert start2 <= end2
		if end1 < start2 or end2 < start1: #Non-overlap
			return 0
		else: 
			return 1
			
			
def find_collisions(sorted_entrez_gene_boundaries):
	"""
	sorted_entrez_gene_boundaries: sorted by start 
	entrez_gene_boundaries:list of (id, start, end)
	return flag_list denoting whether an entry has collisions

	The approach is adopted from merging the islands
	"""
	flag_list = [0] * len(sorted_entrez_gene_boundaries) # entrez id flagged with 1 is to have collision

	current_id = sorted_entrez_gene_boundaries[0][0]
	current_start = sorted_entrez_gene_boundaries[0][1]
	current_end =  sorted_entrez_gene_boundaries[0][2]
	i = 1
	while i < len(sorted_entrez_gene_boundaries):
		compare_id = sorted_entrez_gene_boundaries[i][0]
		compare_start = sorted_entrez_gene_boundaries[i][1]
		compare_end =  sorted_entrez_gene_boundaries[i][2]
		if compare_start > current_end: # not overlapping
			current_id = compare_id
			current_start = compare_start
			current_end = compare_end
		else: #overlapping
			flag_list[i-1] = 1
			flag_list[i] = 1
			current_end = max(current_end, compare_end)
		i += 1
	return flag_list


def test_find_collisions():
	sorted_entrez_gene_boundaries = [(0,1,5), (1,5,6),(2,7,12), (4,20,35), (5,21,70),(6,22,23), (7,100,200)]
	print sorted_entrez_gene_boundaries
	print find_collisions(sorted_entrez_gene_boundaries)


#def annotate_element(my_entrezs, position, chrom):
	#"""
	#my_entrez_collection is an KnownEntrezGenes class
	#element is represented  by a position
	
	#"""
	#my_entrez_dic = my_entrezs.subset_by_chrom(chrom) # return a entrez_genes dic
	#for myid in my_entrez_dic.keys():
		


def remove_entrez_ids_in_collisions (my_entrez_collection, chroms):
	"""
	my_entrez_collection is an KnownEntrezGenes class
	"""
	collison_removed={}
	for chrom in chroms:
		collison_removed.update(remove_entrez_ids_in_collisions_by_chrom(my_entrez_collection, chrom))
	print "There are %d Entrez IDs after removing genomic collisions between entrez_ids." % len(collison_removed.keys())
	# output the entrez ids for those genes in collision
	#all_entrez_ids = get_entrez_ids(entrez_id_collection)
	#remained_entrez_ids = get_entrez_ids(collison_removed_entrez_id_collection)
	#entrez_ids_in_collision=list(set(all_entrez_ids) - set(remained_entrez_ids))
	#output = open('entrez_ids_in_collision', 'w')
	#for myid in entrez_ids_in_collision:
	#	output.write (str(myid) + '\n')
	#output.close()
	return collison_removed
	
def remove_entrez_ids_in_collisions_by_chrom (my_entrez_collection, chrom):
	"""
	my_entrez_collection is a KnownEntrezGenes class
	
	Two entrez IDs could overlap on genome, for these genes, the intron-exon designation could be ambiguous, for analysis of intron retention and 3UTR, it is safer to remove these type of genes 
	
	We assume that within each entrez gene, the refseq genes are already made sure to be overlapping and on the same strand. 
	
	return a dictionary {id:entrez_gene}
	
	"""
	entrez_gene_boundaries=[] # store the id, start, end of each entrez id. so that sortable
	for entrez_id in my_entrez_collection.entrez_ids_by_chrom[chrom]:
		myboundaries = (my_entrez_collection.entrez_genes)[entrez_id].boundaries[0]
		start = myboundaries[0] # the genomic start irrespective of the direction of the gene
		end = myboundaries[1] # the genomic end irrespective of the direction of the gene
		entrez_gene_boundaries.append((entrez_id, start, end))
	sorted_entrez_gene_boundaries = sorted(entrez_gene_boundaries, key = itemgetter(1))	
	flag_list = find_collisions(sorted_entrez_gene_boundaries)
	my_entrez_genes={}
	for i in xrange(len(flag_list)):
		if flag_list[i] == 0: # no collision
			myid = sorted_entrez_gene_boundaries[i][0]
			my_entrez_genes[myid] = (my_entrez_collection.entrez_genes)[myid]
	return my_entrez_genes # return a dictionary keyed by entrez_id
	
def remove_same_strand_entrez_ids_collisions (my_entrez_collection, chroms):
	"""
	my_entrez_collection is a KnownEntrezGenes class
	returns: {id:EntrezGene}
	"""
	collison_removed={}
	for chrom in chroms:
		collison_removed.update(remove_same_strand_entrez_ids_collisions_by_chrom(my_entrez_collection, chrom))
	print "There are ", len(collison_removed.keys()), " Entrez IDs after removing same-strand genomic collisions between entrez_ids"
	return collison_removed

def remove_same_strand_entrez_ids_collisions_by_chrom (my_entrez_collection, chrom):
	"""
	my_entrez_collection is a KnownEntrezGenes class object
	
	Two entrez IDs on the same strand could overlap on genome, for these genes, the intron exon reads designation could be problematic, for analysis of intron retention and 3UTR using strand specific library, it is safer to remove these type of genes 
	
	We assume that within each entrez gene, the refseq genes are already made sure to be overlapping and on the same strand. 
	
	returns: {id:EntrezGene}
	"""
	p_entrez_gene_boundaries=[] # store the id, start, end of each entrez id.
	n_entrez_gene_boundaries=[] 
	for entrez_id in my_entrez_collection.entrez_ids_by_chrom[chrom]:
		# merge the genomic regions defined by refseq ids to find the start and end of the entrez gene 
		myboundaries = (my_entrez_collection.entrez_genes)[entrez_id].boundaries [0]
		start = myboundaries[0] # the genomic start irrespective of the direction of the gene
		end = myboundaries[1] # the genomic end irrespective of the direction of the gene
		refseqs = (my_entrez_collection.entrez_genes)[entrez_id].transcripts
		if refseqs[0].strand == '+':
			p_entrez_gene_boundaries.append([entrez_id, start, end])
		elif refseqs[0].strand == '-':
			n_entrez_gene_boundaries.append([entrez_id, start, end])
	
	# sort entrez_gene_boundaries according to start
	sorted_entrez_gene_boundaries = sorted(p_entrez_gene_boundaries, key = itemgetter(1))	
	flag_list = find_collisions(sorted_entrez_gene_boundaries)
	my_entrez_genes={}
	for i in xrange(len(flag_list)):
		if flag_list[i] == 0:
			myid = sorted_entrez_gene_boundaries[i][0]
			my_entrez_genes[myid] = (my_entrez_collection.entrez_genes)[myid]
	
	# sort entrez_gene_boundaries according to start
	sorted_entrez_gene_boundaries = sorted(n_entrez_gene_boundaries, key = itemgetter(1))		
	flag_list = find_collisions(sorted_entrez_gene_boundaries)
	for i in xrange(len(flag_list)):
		if flag_list[i] == 0:
			myid = sorted_entrez_gene_boundaries[i][0]
			my_entrez_genes[myid] = (my_entrez_collection.entrez_genes)[myid]
	return my_entrez_genes # return a dictionary keyed by entrez_id

def test_gene_structure(entrez_gene_collection, entrez_id):

	#test getExons
	if entrez_id in entrez_gene_collection.entrez_ids:
		gene = entrez_gene_collection.entrez_genes[entrez_id]
		print "Chrom: ", gene.chrom
		print "Entrez ID: ", entrez_id
		print "Strand:", gene.strand
		print "Concordance:", gene.check_concordance()
		print "Boundaries:", gene.boundaries
		print "Number of transcripts:", len(gene.transcripts)
		print
		for transcript in gene.transcripts:
			print
			print "transcript.getAll()"
			print transcript.getAll()
			print 
			print "transcript.getExons()"
			print transcript.getExons()
			print
			print "transcript.getIntrons()"
			print transcript.getIntrons()
			print 
		print 
		print "gene.shared_exonic_regions"
		print gene.shared_exonic_regions
		print "gene.shared_intronic_regions"
		print gene.shared_intronic_regions
	else:
		print entrez_id, " is not among the ids used." 

def test_pickle(pkl):
	annotation = open(pkl, 'rb')
	new_collection = KnownEntrezGenes(chroms, pickle.load(annotation)) 
	annotation.close()
	print new_collection.chroms
	print new_collection.num_genes
	id = new_collection.entrez_genes.keys()[0]
	print id
	gene = new_collection.entrez_genes[id]
	print gene.strand
	for i in new_collection.entrez_genes[id].transcripts:
		print i.getAll()
		
def main(argv):
	parser = OptionParser()
	parser.add_option("-u", "--refseqfile", action="store", type="string", dest="refseqfile", help="input ucsc file for annotated genes, eg,  refFlat_hg19_EntrezID.ucsc", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output the annotation system in pickle", metavar="<file>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);


	
	#Initiate the clustering of transcripts by entrez ID. remove the discordant transcripts (non-overlapping, on opposite strands)
	my_entrez_collection = KnownEntrezGenes.initiate_from_file(chroms, opt.refseqfile)
	#test(my_entrez_collection)
	
	# test module
	test = 1
	if test == 1:
		test_id = my_entrez_collection.entrez_ids[5]
		print "Testing gene structure"
		test_gene_structure(my_entrez_collection, test_id)
		#test_id = 201176
		#if test_id in my_entrez_collection.get_entrez_ids():
			#print "201176 is problematic"
	
	#To ensure contamination-free, first remove collision and then remove NR items because genes may collide with NR

	# Remove all collisions
	print
	print "Remove all collisions"
	collison_removed_entrez_id_collection = KnownEntrezGenes(chroms, remove_entrez_ids_in_collisions(my_entrez_collection, chroms))
	eligible_ids = collison_removed_entrez_id_collection.filter_by_eligibility(collison_removed_entrez_id_collection.entrez_ids, extension)
	collison_removed_entrez_id_collection.output_pickle( opt.outfile + "_collisonremoved.pkl", eligible_ids)
	
	#test_pickle(opt.outfile + "_collisonremoved.pkl")
	#test_id = 201176
	#if test_id in collison_removed_entrez_id_collection.get_entrez_ids():
		#print "201176 is problematic in collison_removed_entrez_id_collection"
	
	
	print 
	print "Remove same-strand collisions"
	#Remove same strand collision
	collison_removed_entrez_id_collection = KnownEntrezGenes(chroms, remove_same_strand_entrez_ids_collisions(my_entrez_collection, chroms))
	eligible_ids = collison_removed_entrez_id_collection.filter_by_eligibility(collison_removed_entrez_id_collection.entrez_ids, extension)
	collison_removed_entrez_id_collection.output_pickle(opt.outfile + "_samestrand_collisonremoved.pkl", eligible_ids)
	#test_pickle(opt.outfile + "_samestrand_collisonremoved.pkl")
	#test_id = 201176
	#if test_id in collison_removed_entrez_id_collection.get_entrez_ids():
		#print "201176 is problematic in same_strand_collison_removed_entrez_id_collection"	

	print 
	print "Remove same-strand collisions and retain only simple genes"
	simple_gene_ids = collison_removed_entrez_id_collection.find_simple_genes(eligible_ids)
	collison_removed_entrez_id_collection.output_pickle(opt.outfile + "_samestrand_collisonremoved_simplegene.pkl", simple_gene_ids)
	#test_pickle(opt.outfile + "_samestrand_collisonremoved_simplegene.pkl")
	
	
	# Examine the histogram of the Entrez vs refseq conversion
	collison_removed_entrez_id_collection.id_clustering_histogram("Entrezid_refseq_NM_NR_ids.hist")
	
	
	
if __name__ == "__main__":
	main(sys.argv)
