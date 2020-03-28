import HTSeq
import sys, os.path, re
from string import *
from optparse import OptionParser
import collections

def find_TSS_region_iv(transcript_iv, TSS_region_length):
        TSS_region_iv = HTSeq.GenomicInterval(transcript_iv.chrom, transcript_iv.start_d - TSS_region_length/2, transcript_iv.start_d + TSS_region_length/2, transcript_iv.strand)
        return TSS_region_iv

def find_promoter_iv(transcript_iv, promoter_upstream_extension, promoter_downstream_extension):
        if transcript_iv.strand == "+":
                promoter_iv = HTSeq.GenomicInterval(transcript_iv.chrom, transcript_iv.start_d - promoter_upstream_extension, transcript_iv.start_d + promoter_downstream_extension, "+")
        elif transcript_iv.strand == "-":
                promoter_iv = HTSeq.GenomicInterval(transcript_iv.chrom, transcript_iv.start_d - promoter_downstream_extension, transcript_iv.start_d + promoter_upstream_extension, "-")
        return promoter_iv

def main(argv):
        parser = OptionParser()
        parser.add_option("-i", "--input_bed_file", action="store", type="string", dest="inputfile", metavar="<file>", help="input bed file that records peaks")
        parser.add_option("-g", "--genes_GTF_file", action="store", type="string", dest="genesfile", metavar="<file>", help="GTF genes annotation file")
        parser.add_option("-u", "--promoter_upstream_extension", action="store", type="int", dest="promoter_upstream_extension", metavar="<int>", help="promoter upstream extension")
        parser.add_option("-d", "--promoter_downstream_extension", action="store", type="int", dest="promoter_downstream_extension", metavar="<int>", help="promoter downstream extension")   
        parser.add_option("-t", "--TSS_region_length", action="store", type="int", dest="TSS_region_length", metavar="<int>", help="TSS region length")
        parser.add_option("-o", "--output_distribution_file", action="store", type="string", dest="outfile", metavar="<file>", help="output peaks genic region distribution file")     
        (opt, args) = parser.parse_args(argv)
        if len(argv) < 10:
                parser.print_help()
                sys.exit(1)

        promoter_upstream_extension = opt.promoter_upstream_extension    
        if not opt.promoter_downstream_extension:
                promoter_downstream_extension = - opt.TSS_region_length / 2
        else:
                promoter_downstream_extension = opt.promoter_downstream_extension
        TSS_region_length = opt.TSS_region_length


        TSS_region = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        promoters = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        gene_body = HTSeq.GenomicArrayOfSets("auto", stranded=False) 

        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)

        for feature in gtffile:
                if feature.type == "transcript_region":
                        TSS_region_iv = find_TSS_region_iv(feature.iv, TSS_region_length)
                        TSS_region[TSS_region_iv] += ("TSS_region", feature.attr["gene_id"])
                        promoter_iv = find_promoter_iv(feature.iv, promoter_upstream_extension, promoter_downstream_extension)
                        promoters[promoter_iv] += ("promoter", feature.attr["gene_id"])
                elif feature.type == "gene_region":
                        gene_body[feature.iv] += ("gene_body", feature.attr["gene_id"])    
                        
        peak_genic_region = collections.defaultdict(lambda: set())		
        bedfile = HTSeq.BED_Reader(opt.inputfile)   
        
        TSS_region_peaks_set = set()
        promoter_peaks_set = set()
        gene_body_peaks_set = set()
        intergenic_region_peaks_set = set()
        
        for peak in bedfile:
                for iv, step_set in TSS_region[peak.iv].steps():
                        if step_set:
                                TSS_region_peaks_set.add(peak.iv)
                                for region_type in list(step_set):
                                        peak_genic_region[peak.iv].add(region_type)
                if peak.iv not in peak_genic_region:
                        peak_genic_region[peak.iv] = set()     
                        
        for peak_iv in peak_genic_region:
                if len(peak_genic_region[peak_iv]) == 0:
                        for iv,step_set in promoters[peak_iv].steps():
                                if step_set:
                                        promoter_peaks_set.add(peak_iv)
                                        for region_type in list(step_set):
                                                peak_genic_region[peak_iv].add(region_type)       
                                                
        for peak_iv in peak_genic_region:
                if len(peak_genic_region[peak_iv]) == 0:
                        for iv,step_set in gene_body[peak_iv].steps():
                                if step_set:
                                        gene_body_peaks_set.add(peak_iv)
                                        for region_type in list(step_set):
                                                peak_genic_region[peak_iv].add(region_type)     
                                                
        for peak_iv in peak_genic_region:
                if len(peak_genic_region[peak_iv]) == 0:
                        intergenic_region_peaks_set.add(peak_iv)
                        peak_genic_region[peak_iv].add("intergenic_region")     
                        
        outfile = open(opt.outfile, "w")
        header = "#" + "chrom" + "\t" + "start" + "\t" + "end" + "\t" + "genic_region" + "\n"
        outfile.write(header)
        for peak_iv in sorted(peak_genic_region.keys(), key=lambda f: ( f.chrom, f.start )):
                if len(peak_genic_region[peak_iv]) > 1:
                        region_type_list = []
                        for region_type in list(peak_genic_region[peak_iv]):
                                region_type_list.append(" ".join(region_type))
                        region_type_list_str = ",".join(region_type_list)

                        outline = peak_iv.chrom + "\t" + str(peak_iv.start) + "\t" + str(peak_iv.end) + "\t" + region_type_list_str + "\n"
                        outfile.write(outline)

                elif len(peak_genic_region[peak_iv]) == 1:
                        if len(list(peak_genic_region[peak_iv])[0]) == 2:
                                outline = peak_iv.chrom + "\t" + str(peak_iv.start) + "\t" + str(peak_iv.end) + "\t" + " ".join(list(peak_genic_region[peak_iv])[0]) + "\n"
                                outfile.write(outline)
                        elif len(list(peak_genic_region[peak_iv])[0]) == 1:
                                outline = peak_iv.chrom + "\t" + str(peak_iv.start) + "\t" + str(peak_iv.end) + "\t" + list(peak_genic_region[peak_iv])[0] + "\n"
                                outfile.write(outline)	

        print "Number of peaks in TSS region is %d, number of peaks in promoter is %d, number of peaks in gene body is %d, and number of peaks in intergenic region is %d." % (len(TSS_region_peaks_set), len(promoter_peaks_set), len(gene_body_peaks_set), len(intergenic_region_peaks_set))

if __name__ == "__main__":
        main(sys.argv)

