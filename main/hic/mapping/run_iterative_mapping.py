import sys
from optparse import OptionParser
from hiclib.mapping import iterative_mapping
import HTSeq

def main(argv):
        parser = OptionParser()	
        parser.add_option("-b", "--bowtie_path", action="store", type="string", dest="bowtie_path", metavar="<file>", help="bowtie path")
        parser.add_option("-i", "--bowtie_index_path", action="store", type="string", dest="bowtie_index_path", metavar="<file>", help="bowtie index path")
	parser.add_option("-q", "--fastq_path", action="store", type="string", dest="fastq_path", metavar="<file>", help="fastq path")
	parser.add_option("-o", "--out_sam_path", action="store", type="string", dest="out_sam_path", metavar="<file>", help="output sam path")
	parser.add_option("-m", "--min_seq_len", action="store", type="int", dest="min_seq_len", metavar="<int>", help="The truncation length at the first iteration of mapping")
	parser.add_option("-t", "--len_step", action="store", type="int", dest="len_step", metavar="<int>", help="The increase in truncation length at each iteration")
	parser.add_option("-s", "--seq_start", action="store", type="int", dest="seq_start", metavar="<int>", help="seq start")
	parser.add_option("-e", "--seq_end", action="store", type="int", dest="seq_end", metavar="<int>", help="seq end")
        (opt, args) = parser.parse_args(argv)
        if len(argv) < 12:
                parser.print_help()
                sys.exit(1)
	
	if opt.seq_start and opt.seq_end:
		#iterative_mapping(opt.bowtie_path, opt.bowtie_index_path, opt.fastq_path, opt.out_sam_path, opt.min_seq_len, opt.len_step, seq_start=opt.seq_start, seq_end=opt.seq_end)
		iterative_mapping(opt.bowtie_path, opt.bowtie_index_path, opt.fastq_path, opt.out_sam_path, opt.min_seq_len, opt.len_step, seq_start=opt.seq_start, seq_end=opt.seq_end, bowtie_flags='--no-unal --mm')
	else:
		#iterative_mapping(opt.bowtie_path, opt.bowtie_index_path, opt.fastq_path, opt.out_sam_path, opt.min_seq_len, opt.len_step)
		iterative_mapping(opt.bowtie_path, opt.bowtie_index_path, opt.fastq_path, opt.out_sam_path, opt.min_seq_len, opt.len_step, bowtie_flags='--no-unal --mm')
        
if __name__ == "__main__":
	main(sys.argv)
                        
