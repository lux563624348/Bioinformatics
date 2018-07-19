### XIANG LI

import HTSeq
import sys, os
from optparse import OptionParser
import numpy
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
## FUNCTIONS
def Get_Site_Profile(read_file, gtf_file, gene_list, site_name, window_size, fragment_size, resolution, upstreamExtension,downstreamExtension):
    ## Read Reads_File, and distribute them into GenomicArray.
    num_reads = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    bedfile = HTSeq.BED_Reader(read_file)
    for alt in bedfile:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        ga[alt_pos] += 1
        num_reads += 1
##########################################################################################
    ## Read gtf_file
    gtffile = HTSeq.GFF_Reader(gtf_file)
##########################################################################################    
    ## Read Gene_list_file
    gene_list = pd.read_csv(gene_list, header=None, sep='\t', usecols=[3])[3].str.strip()
    gene_list_unique=gene_list.unique()
##########################################################################################
    ### Set up number of points 
    upstream_num_points = upstreamExtension / resolution
    downstream_num_points = downstreamExtension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    ### To calculate the number of transcripts that for profile plot.
    ### And also store position in a set.
    num_transcripts=0
    site_pos_set = set()
    
    if (site_name=='TSS'):
        for feature in gtffile:
            if feature.type == "exon" and feature.attr["gene_id"] in gene_list_unique:
                site_pos_set.add(feature.iv.start_d_as_pos)
                num_transcripts += 1
        print num_transcripts
    elif(site_name=='TES'):
        for feature in gtffile:
            if feature.type == "exon" and feature.attr["gene_id"] in gene_list_unique:
                site_pos_set.add(feature.iv.end_d_as_pos)
                num_transcripts += 1
        print num_transcripts
    else:
        print("ERROR: Default is TSS or TES, your input is neither!")
    
    
    profile = numpy.zeros(total_num_points, numpy.int)
    
    for tss_pos in site_pos_set:
            index = 0
            while index < total_num_points:
                count_in_window = 0
                if tss_pos.strand == "+":
                    index_pos = tss_pos.pos + (index - upstream_num_points) * resolution
                    tss_pos_window_iv = HTSeq.GenomicInterval(tss_pos.chrom, index_pos - window_size / 2, index_pos + window_size / 2)
                elif tss_pos.strand == "-":
                    index_pos = tss_pos.pos - (index - upstream_num_points) * resolution
                    tss_pos_window_iv = HTSeq.GenomicInterval(tss_pos.chrom, index_pos - window_size / 2 + 1, index_pos + window_size / 2 + 1)

                for step_iv, step_count in ga[tss_pos_window_iv].steps():
                    count_in_window += step_count * step_iv.length
                profile[index] += count_in_window
                index += 1
                
    return (num_reads, num_transcripts, profile)

def profile_normalization(profile, num_reads, num_transcripts, window_size, adjusted_normalization):
    normalization = 1.0
    normalization = num_reads/1000000.0
    normalization *= num_transcripts
    normalization *= window_size/1000.0
    normalization *= adjusted_normalization

    print "Number of locations: %i" % num_transcripts
    print "Number of reads: %i" % num_reads
    print "Normalization is by total number of reads per million. normalization = %f" % normalization	

    return profile / normalization

def profile_plot_site(norm_profile, resolution, upstreamExtension, downstreamExtension, genes_set_name, con_name, site_name):
    ### TSS or TES
    
    fig, ax = plt.subplots(1,1)
    ax.plot(numpy.arange(-upstreamExtension/resolution, downstreamExtension/resolution+1), norm_profile, label=con_name)
    plt.title(genes_set_name)
    plt.legend(loc='upper right')
    
    x=[-upstreamExtension/resolution, 0, downstreamExtension/resolution]
    ax.set_xticks(x)
    customized_xticks=['-'+ str(upstreamExtension),site_name,str(downstreamExtension)]
    ax.set_xticklabels(customized_xticks, fontsize=18)
    ax.grid(which='major', axis='x', linestyle='--')
    fig.savefig('Profile_'+site_name+'_'+con_name+'.png')
    
    ### Also Input Profile txt
    outfile='Profile_'+site_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    xValues = numpy.arange(-upstreamExtension / resolution, downstreamExtension / resolution + 1, 1)
    xValues *= resolution
    for index in range(len(xValues)):
        outline = str(xValues[index]) + "\t" + str(norm_profile[index]) + "\n"
        f.write(outline)
    f.close()

### FUNCTION


### FUNCTIONS
def main(argv):
	desc="""This is a template for the analysis of aggretated tag distribution with respect to a set of points, such as the TSSs of known genes, with one profile from each strand."""
	parser = OptionParser(description=desc)
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bed_file", help="file with tags in bed format", metavar="<file>")
	parser.add_option("-c", "--condition_name", action="store", type="string",
		dest="condition", help="Condition of Input_Reads", metavar="<str>")
	parser.add_option("-g", "--genes_gtf_file", action="store", type="string",
			dest="gtf_file", help="genes gtf file", metavar="<file>")
	parser.add_option("-k", "--known_gene_list", action="store", type="string",
			dest="known_gene_list", help="known gene list", metavar="<file>")
	parser.add_option("-l", "--list_name", action="store", type="string",
		dest="list_name", help="Name of your input gene list", metavar="<str>")
	parser.add_option("-t", "--TypeOfSites", action="store", type="string",
			dest="type_site", help="TSS, TES, TFBS", metavar="<str>")
	parser.add_option("-o", "--outfolder", action="store", type="string",
			dest="outfolder", help="output folder", metavar="<file>")
	parser.add_option("-n", "--normalization", action="store", type="float", dest="norm",
			help="additional normalization in addition to number of sites, number of reads per million and window_size per 1K")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size",
	                  help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps", metavar="<int>")
	parser.add_option("-u", "--UpstreamExtension", action="store", type="int",
			dest="upstreamExtension", help="UpstreamExtension", metavar="<int>")
	parser.add_option("-d", "--DownstreamExtension", action="store", type="int",
			dest="downstreamExtension", help="DownstreamExtension", metavar="<int>")
	parser.add_option("-w", "--WindowSize", action="store", type="int",
			dest="window_size", help="window size for averaging. When window size > resolution, there is smoothing", metavar="<int>")	
	parser.add_option("-r", "--resolution", action="store", type="int",
			dest="resolution", help="resolution of the upstream and downstream profile, eg, 5", metavar="<int>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 26:
		parser.print_help()
		sys.exit(1)
		
	print "Here is the Summary of your input."
	print "Name of Input Read Condition: %s" % opt.condition
	print "Name of Input Gene list: %s" % opt.list_name
	print "Generating Profiles around: %s" % opt.type_site
	print "Upstream extension: %i" % opt.upstreamExtension
	print "Downstream extension: %i" % opt.downstreamExtension
	print "Input Fragment Size: %i" % opt.fragment_size
	print "Scanning window size: %i" % opt.window_size
	print "Scanning resolution: %i" % opt.resolution
	print "Adjusted Normalization: %f" % opt.norm
	print "End of Summary."
	print " "
	
	### Three steps 
	#### First
	(num_reads, num_transcripts, Site_profile) = Get_Site_Profile(opt.bed_file, opt.gtf_file, opt.known_gene_list, opt.type_site, 
	opt.fragment_size, opt.window_size, opt.resolution, opt.upstreamExtension, opt.downstreamExtension)
	#### Second
	norm_Site_profile = profile_normalization (Site_profile, num_reads, num_transcripts, opt.window_size, opt.norm)
	#### Third
	profile_plot_site(norm_Site_profile, opt.resolution, opt.upstreamExtension, opt.downstreamExtension, opt.list_name, opt.condition, opt.type_site)

if __name__ == "__main__":
	main(sys.argv)

