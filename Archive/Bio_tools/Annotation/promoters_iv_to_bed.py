import HTSeq
import sys, os, re
from string import *
from optparse import OptionParser
import collections
import itertools

def main(argv):
    parser = OptionParser()
    parser.add_option("-g", "--genes_GTF_file", action="store", type="string", dest="genesfile", metavar="<file>", help="GTF genes annotation file")
    parser.add_option("-k", "--gene_list", action="store", type="string", dest="genelist", metavar="<file>", help="gene list")
    parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'PromoterGenebody' or 'NonPromoter' or 'TSS' or 'three_UTR_region'", action="store", type="string", dest="region_type", metavar="<str>", help="region")
    parser.add_option("-u", "--promoter_upstream_extension", action="store", type="int", dest="promoter_upstream_extension", help="upstream extension of promoter region from TSS", metavar="<int>")
    parser.add_option("-d", "--promoter_downstream_extension", action="store", type="int", dest="promoter_downstream_extension", help="downstream extension of promoter region from TSS",  metavar="<int>")
    parser.add_option("-s", "--stranded", action="store_true", dest="stranded", metavar="<file>", help="output stranded interval")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name containg all region intervals")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)    
        
    if opt.promoter_upstream_extension or opt.promoter_upstream_extension == 0:
        UpstreamExtension = opt.promoter_upstream_extension
    if opt.promoter_downstream_extension or opt.promoter_downstream_extension == 0:
        DownstreamExtension = opt.promoter_downstream_extension   
    
    outfile = open(opt.outfile, 'w')    
    
    if opt.region_type == "Promoter":
        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
        if opt.genelist:
            gene_list = []
            with open(opt.genelist, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    gene_list.append(gene_id)
            gtffile = filter(lambda feature: feature.attr['gene_id'] in gene_list, gtffile)
            
        for feature in gtffile:
            if feature.type == "transcript_region":
                transcript_region_iv = feature.iv
                if transcript_region_iv.strand == "+":                   
                    if transcript_region_iv.start_d - UpstreamExtension < 0:
                        promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.start_d + DownstreamExtension + 1, transcript_region_iv.strand)
                    else:
                        promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start_d - UpstreamExtension, transcript_region_iv.start_d + DownstreamExtension + 1, transcript_region_iv.strand)
                elif transcript_region_iv.strand == "-":
                    if transcript_region_iv.start_d - DownstreamExtension < 0:
                        promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.start_d + UpstreamExtension + 1, transcript_region_iv.strand)
                    else:
                        promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start_d - DownstreamExtension, transcript_region_iv.start_d + UpstreamExtension + 1, transcript_region_iv.strand)
                
                if opt.stranded:
                    outline = promoter_iv.chrom + '\t' + str(promoter_iv.start) + '\t' + str(promoter_iv.end) + '\t' + feature.attr['gene_id'] + '\t' + '0' + '\t' + promoter_iv.strand + '\n'
                else:
                    outline = promoter_iv.chrom + '\t' + str(promoter_iv.start) + '\t' + str(promoter_iv.end) + '\t' + feature.attr['gene_id'] + '\n'
                outfile.write(outline)
                
    elif opt.region_type == "PromoterGenebody":
        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
        if opt.genelist:
            gene_list = []
            with open(opt.genelist, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    gene_list.append(gene_id)
            gtffile = filter(lambda feature: feature.attr['gene_id'] in gene_list, gtffile)
            
        for feature in gtffile:
            if feature.type == "transcript_region":  
                transcript_region_iv = feature.iv
                if transcript_region_iv.strand == "+":
                    if transcript_region_iv.start - UpstreamExtension < 0:
                        promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.end + DownstreamExtension, transcript_region_iv.strand)
                    else:
                        promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start - UpstreamExtension, transcript_region_iv.end + DownstreamExtension, transcript_region_iv.strand)
                elif transcript_region_iv.strand == "-":
                    if transcript_region_iv.start_d - DownstreamExtension < 0:
                        promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0, transcript_region_iv.end + UpstreamExtension, transcript_region_iv.strand)
                    else:    
                        promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start - DownstreamExtension, transcript_region_iv.end + UpstreamExtension, transcript_region_iv.strand)
                
                if opt.stranded:
                    outline = promoter_plus_genebody_iv.chrom + '\t' + str(promoter_plus_genebody_iv.start) + '\t' + str(promoter_plus_genebody_iv.end) + '\t' + feature.attr['gene_id'] + '\t' + '0' + '\t' + promoter_plus_genebody_iv.strand + '\n'
                else:
                    outline = promoter_plus_genebody_iv.chrom + '\t' + str(promoter_plus_genebody_iv.start) + '\t' + str(promoter_plus_genebody_iv.end) + '\t' + feature.attr['gene_id'] + '\n'
                outfile.write(outline)
    
    elif opt.region_type == "NonPromoter":
        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
        if opt.genelist:
            gene_list = []
            with open(opt.genelist, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    gene_list.append(gene_id)
            gtffile = filter(lambda feature: feature.attr['gene_id'] in gene_list, gtffile)
        
        for feature in gtffile:
            if feature.type == "transcript_region":  
                transcript_region_iv = feature.iv
                if transcript_region_iv.length <= UpstreamExtension:
                    continue
                if transcript_region_iv.strand == "+":
                    non_promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start + UpstreamExtension, transcript_region_iv.end + DownstreamExtension, transcript_region_iv.strand)
                elif transcript_region_iv.strand == "-":
                    non_promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start - DownstreamExtension, transcript_region_iv.end - UpstreamExtension, transcript_region_iv.strand)
                
                if opt.stranded:
                    outline = non_promoter_iv.chrom + '\t' + str(non_promoter_iv.start) + '\t' + str(non_promoter_iv.end) + '\t' + feature.attr['gene_id'] + '\t' + '0' + '\t' + non_promoter_iv.strand + '\n'
                else:
                    outline = non_promoter_iv.chrom + '\t' + str(non_promoter_iv.start) + '\t' + str(non_promoter_iv.end) + '\t' + feature.attr['gene_id'] + '\n'
                outfile.write(outline)        
        
    elif opt.region_type == "TSS":
        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
        if opt.genelist:
            gene_list = []
            with open(opt.genelist, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    gene_list.append(gene_id)
            gtffile = filter(lambda feature: feature.attr['gene_id'] in gene_list, gtffile)  


        for feature in gtffile:
            if feature.type == "transcript_region":  
                transcript_region_iv = feature.iv
                TSS_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, transcript_region_iv.start_d, transcript_region_iv.start_d + 1, transcript_region_iv.strand)
 
                if opt.stranded:
                    outline = TSS_iv.chrom + '\t' + str(TSS_iv.start) + '\t' + str(TSS_iv.end) + '\t' + feature.attr['gene_id'] + '\t' + '0' + '\t' + TSS_iv.strand + '\n'
                else:
                    outline = TSS_iv.chrom + '\t' + str(TSS_iv.start) + '\t' + str(TSS_iv.end) + '\t' + feature.attr['gene_id'] + '\n'
                outfile.write(outline)
                
    elif opt.region_type == "three_UTR_region":
        gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
        if opt.genelist:
            gene_list = []
            with open(opt.genelist, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    gene_list.append(gene_id)
            gtffile = filter(lambda feature: feature.attr['gene_id'] in gene_list, gtffile)
            
        for feature in gtffile:
            if feature.type == "three_UTR_region":  
                three_UTR_region_iv = feature.iv
 
                if opt.stranded:
                    outline = three_UTR_region_iv.chrom + '\t' + str(three_UTR_region_iv.start) + '\t' + str(three_UTR_region_iv.end) + '\t' + feature.attr['gene_id'] + '\t' + '0' + '\t' + three_UTR_region_iv.strand + '\n'
                else:
                    outline = three_UTR_region_iv.chrom + '\t' + str(three_UTR_region_iv.start) + '\t' + str(three_UTR_region_iv.end) + '\t' + feature.attr['gene_id'] + '\n'
                outfile.write(outline)        
            
                
if __name__ == "__main__":
    main(sys.argv)
            
            
