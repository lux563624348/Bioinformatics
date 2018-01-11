import re, os, sys, shutil
from string import *
from optparse import OptionParser
import GenomeData

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	(opt, agrs) = parser.parse_args(argv)
	
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	if opt.species in GenomeData.species_chroms.keys():
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
		max = 0
		maxchrom = ''
		for chrom in GenomeData.species_chroms[opt.species]:
			if chrom_lengths[chrom] > max:
				max = chrom_lengths[chrom]
				maxchrom = chrom
		print opt.species, maxchrom, max
	else:
		print opt.species + " is not recognized"
		

if __name__ == "__main__":
	main(sys.argv)




