########################################################################
## 07/20/2018
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
########################################################################
# Usage 
#python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} \
#	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
########################################################################

import HTSeq
import sys, re
from optparse import OptionParser
import numpy
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
## FUNCTIONS
### FUNCTION
def Get_Site_Profile(INPUT_read_file, INPUT_gtf_file, INPUT_gene_list, site_name, window_size, fragment_size, resolution, upstreamExtension,downstreamExtension):
    ## Read Reads_File, and distribute them into GenomicArray.
    num_reads = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    bedfile = HTSeq.BED_Reader(INPUT_read_file)
    for alt in bedfile:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        ga[alt_pos] += 1
        num_reads += 1
##########################################################################################
    ## Read gtf_file
    gtffile = HTSeq.GFF_Reader(INPUT_gtf_file)
##########################################################################################    
    ## Read Gene_list_file
    gene_list = pd.read_csv(INPUT_gene_list, header=None, sep='\s+', usecols=[3])[3].str.strip()
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
        #print num_transcripts
    elif(site_name=='TES'):
        for feature in gtffile:
            if feature.type == "exon" and feature.attr["gene_id"] in gene_list_unique:
                site_pos_set.add(feature.iv.end_d_as_pos)
                num_transcripts += 1
        #print num_transcripts
    else:
        print("ERROR: Default is TSS or TES, your input is neither!")
    
    
    profile = numpy.zeros(total_num_points)
    
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
                ###gives average read count for position
                profile[index] += count_in_window / (1.0*window_size)
                
                index += 1
    
    print "...................................................."
    print "Profile on: %s " % site_name
    print "Number of locations: %i" % num_transcripts
    print "Number of reads: %i" % num_reads 
            
    return  ((profile*(10**9))/(num_reads*num_transcripts))

def Get_GeneBody_Profile(INPUT_read_file, INPUT_gtf_file, INPUT_gene_list, genic_partition, fragment_size):
    ###
    bed_file = HTSeq.BED_Reader(INPUT_read_file)
##########################################################################################
    ### rename gene list format and gtf file.
    gene_list = pd.read_csv(INPUT_gene_list, header=None, sep='\s+')
    gene_list = gene_list.rename(columns={0:'chr', 1:'start',2:'end', 3:'gene_id'})
    gene_list['gene_id'] = gene_list['gene_id'].str.strip()


    gtf_file = pd.read_csv(INPUT_gtf_file, header=None, sep='\s+', usecols=[0,3,4,11,6])
    gtf_file = gtf_file.rename(columns={0:'chr', 3:'start', 4:'end', 11:'gene_id', 6:'strand'})
    gtf_file['gene_id']=gtf_file['gene_id'].str.replace(";","").str.strip()
##########################################################################################
    ### Then across reference, merge gtf and genelist, basically, get strand information from gtf, and merged into gene list.
    gene_list_with_strand = pd.merge(gene_list, gtf_file, how='inner', on=['chr', 'start', 'end', 'gene_id'])
    ### A strong merge condition, requires that chr, start, end and gene_id has to be exactly same.
########################################################################################## 
    if (len(gene_list_with_strand) != len(gene_list)):
        print "Input Genelist is not constant with reference GTF, please doublecheck!"
    else:
        print("Looks good so far!_Comment_From_XIANG_LI")
##########################################################################################


## Read Reads_File, and distribute them into GenomicArray.
    num_reads = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")

    for alt in bed_file:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        ga[alt_pos] += 1
        num_reads += 1
##########################################################################################
### After distribution, now we can calculate profile.
    profile = numpy.zeros(genic_partition)

    Num_Skip = 0
    for i in range(len(gene_list_with_strand)):
        partition_size = (gene_list_with_strand.loc[i,'end']-gene_list_with_strand.loc[i,'start']) / (1.0*genic_partition) ## Prevent int division
        if(partition_size < 1):
            Num_Skip +=1 
            continue
##########################################################################################

        site_iv = HTSeq.GenomicInterval(gene_list_with_strand.loc[i,'chr'], gene_list_with_strand.loc[i,'start'],
                                        gene_list_with_strand.loc[i,'end'], gene_list_with_strand.loc[i,'strand'])
        index = 0

        for site_pos in site_iv.xrange_d(partition_size):
            count_in_window = 0
            if site_pos.strand == "+":
                site_pos_window_iv = HTSeq.GenomicInterval(site_pos.chrom, site_pos.pos, site_pos.pos + partition_size)
            elif site_pos.strand == "-":
                site_pos_window_iv = HTSeq.GenomicInterval(site_pos.chrom, site_pos.pos - partition_size + 1,
                                                           site_pos.pos + 1)

            for step_iv, step_count in ga[site_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[index] += count_in_window / (1.0*partition_size)
            index += 1
            if index >= genic_partition:
                break
	print "...................................................."
    print ("Total Number of Skipped gene_list: "+str(Num_Skip))
    print ("Because their length are less than your choose genic_partition")
    plt.plot((profile*(10**9))/(num_reads*(len(gene_list_with_strand)-Num_Skip)))
    plt.show()
    return ( (profile*(10**9))/(num_reads*(len(gene_list_with_strand)-Num_Skip)) )

def profile_plot_site(norm_profile, resolution, upstreamExtension, downstreamExtension, genes_set_name, con_name, site_name):
    ### TSS or TES
    
    fig, ax=plt.subplots(1,1)
    ax.plot(numpy.arange(-upstreamExtension/resolution, downstreamExtension/resolution+1), norm_profile, label=con_name)
    plt.title(genes_set_name)
    plt.legend(loc='upper right')
    
    x=[-upstreamExtension/resolution, 0, downstreamExtension/resolution]
    ax.set_xticks(x)
    customized_xticks=['-'+ str(upstreamExtension),site_name,str(downstreamExtension)]
    ax.set_xticklabels(customized_xticks, fontsize=18)
    ax.grid(which='major', axis='x', linestyle='--')
    plt.show()
    fig.savefig('Profile_'+site_name+'_'+genes_set_name+'_'+con_name+'.png')
    
    ### also Input Profile txt
    outfile='Profile_'+site_name+'_'+genes_set_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    xValues = numpy.arange(-upstreamExtension / resolution, downstreamExtension / resolution + 1, 1)
    for index in range(len(xValues)):
        outline = str(xValues[index]) + "\t" + str(norm_profile[index]) + "\n"
        f.write(outline)
    f.close()
    return None

def profile_Up_genebody_Down_site(upstream_profile, site_body_profile, downstream_profile, resolution, genic_partition,
                                  upstreamExtension, downstreamExtension, genes_set_name, con_name, site_name):
    ### Upstream Profile Saving
    outfile='Profile_TSS_'+genes_set_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    header = "#upstream_region" + "\n"
    f.write(header)
    upstream_xValues = numpy.arange(0, upstreamExtension + 1, resolution)[-1::-1] * (-1)
    for index in range(len(upstream_xValues)):
        outline = str(upstream_xValues[index]) + "\t" + str(upstream_profile[index]) + "\n"
        f.write(outline)
    f.close()
##########################################################################################
    ### Downstream Profile Saving
    outfile='Profile_TES_'+genes_set_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    header = "#downstream_region" + "\n"
    f.write(header)
    downstream_xValues = numpy.arange(0, downstreamExtension + 1, resolution)
    for index in range(len(downstream_xValues)):
        outline = str(downstream_xValues[index]) + "\t" + str(downstream_profile[index]) + "\n"
        f.write(outline)
    f.close()
##########################################################################################
    ### Genebody Profile Saving
    outfile='Profile_SiteBody_'+genes_set_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    header = "#SiteBody_region" + "\n"
    f.write(header)
    site_body_xValues = [0.0] * genic_partition
    for i in xrange(genic_partition):
        site_body_xValues[i] = (i + 0.5) / genic_partition
    for index in range(len(site_body_xValues)):
        outline = str(site_body_xValues[index]) + "\t" + str(site_body_profile[index]) + "\n"
        f.write(outline)  
    f.close()
##########################################################################################    
    ### Plot     
    numTicksInBody = 5

    upstream_scale = 0.5
    upstream_scaled_xValues = [numTicksInBody * upstream_scale * (1 + item * 1.0 / upstreamExtension) for item in
                               upstream_xValues]

    site_body_scaled_xValues = [numTicksInBody * upstream_scale + item * numTicksInBody for item in
                                site_body_xValues]

    downstream_scale = 0.5
    downstream_scaled_xValues = [
        numTicksInBody * upstream_scale + numTicksInBody + numTicksInBody * downstream_scale * item * 1.0 / downstreamExtension
        for item in downstream_xValues]

    xValues = upstream_scaled_xValues + site_body_scaled_xValues + downstream_scaled_xValues
    
    profile = numpy.append(upstream_profile[0:upstreamExtension/resolution+1], site_body_profile)
    profile = numpy.append(profile, downstream_profile[downstreamExtension/resolution:])

    xticks_subset = [0] * (numTicksInBody + 1)
    xticklabels_subset = [0] * (numTicksInBody + 1)
    for i in xrange(numTicksInBody + 1):
        xticks_subset[i] = numTicksInBody * upstream_scale + i
        if i == 0:
            xticklabels_subset[i] = 'Start'
        elif i == numTicksInBody:
            xticklabels_subset[i] = 'End'
        else:
            xticklabels_subset[i] = str(i * 1.0 / numTicksInBody)

    xticks = [0] + xticks_subset + [xticks_subset[-1] + numTicksInBody * downstream_scale]
    xticklabels = [ '-'+str(upstreamExtension)] + xticklabels_subset + [
        str(downstreamExtension)]


    fig, ax=plt.subplots(1,1, figsize=(12,6))
    ax.plot(xValues, profile, "b", label=con_name)
    ax.grid(which='major', axis='x', linestyle='--')
    ax.set_title(genes_set_name)
    plt.xlabel('Site Coordinate', fontsize=12)
    plt.legend(loc='upper right')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=14)
    plt.ylabel('RPKM', fontsize=24)
    fig_name='Profile-SiteBody-'+genes_set_name+'_'+con_name
    title = re.sub("_", " ", fig_name)
    plt.title(title, fontsize=16)
    #plt.savefig(fig_name + '_plot.eps', format='eps')
    plt.savefig(fig_name + '.png', format='png')

    return None
    ### TSS or TES
    
    fig, ax=plt.subplots(1,1)
    ax.plot(numpy.arange(-upstreamExtension/resolution, downstreamExtension/resolution+1), norm_profile, label=con_name)
    plt.title(genes_set_name)
    plt.legend(loc='upper right')
    
    x=[-upstreamExtension/resolution, 0, downstreamExtension/resolution]
    ax.set_xticks(x)
    customized_xticks=['-'+ str(upstreamExtension),site_name,str(downstreamExtension)]
    ax.set_xticklabels(customized_xticks, fontsize=18)
    ax.grid(which='major', axis='x', linestyle='--')
    plt.show()
    fig.savefig('Profile_'+site_name+'_'+con_name+'.png')
    
    ### also Input Profile txt
    outfile='Profile_'+site_name+'_'+con_name+'.txt'
    f = open(outfile, "w")
    xValues = numpy.arange(-upstreamExtension / resolution, downstreamExtension / resolution + 1, 1)
    for index in range(len(xValues)):
        outline = str(xValues[index]) + "\t" + str(norm_profile[index]) + "\n"
        f.write(outline)
    f.close()


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
			dest="type_site", help="TSS, TES, TFBS, GeneBody", metavar="<str>")
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
	parser.add_option("-p", "--partition_points", action="store", type="int",
			dest="genic_partition", help="genic_partition of genebody, eg, 10", metavar="<int>", default=10)
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 28:
		parser.print_help()
		sys.exit(1)
	
	print " "
	print "Here is the Summary of your input."
	print "Name of Input Read Condition: %s" % opt.condition
	print "Name of Input Gene list: %s" % opt.list_name
	print "Generating Profiles around: %s" % opt.type_site
	print "Genebody Genic_Partition Points(Default: 10) %i" genic_partition
	print "Upstream extension: %i" % opt.upstreamExtension
	print "Downstream extension: %i" % opt.downstreamExtension
	print "Input Fragment Size: %i" % opt.fragment_size
	print "Scanning window size: %i" % opt.window_size
	print "Scanning resolution: %i" % opt.resolution
	print "Adjusted Normalization: %f" % opt.norm
	print "End of Summary."
	print " "
	
	### Two steps
	
	#### First GeneBoydy
	if (opt.type_site == 'GeneBody' ):
		###Up
		Upprofile = Get_Site_Profile(opt.bed_file, opt.gtf_file, opt.known_gene_list, 'TSS', opt.window_size,
		opt.fragment_size, opt.resolution, opt.upstreamExtension, opt.downstreamExtension)
		###Down
		Downprofile = Get_Site_Profile(opt.bed_file, opt.gtf_file, opt.known_gene_list, 'TSS', opt.window_size,
		opt.fragment_size, opt.resolution, opt.upstreamExtension, opt.downstreamExtension)
		### GENEBODY
		genebody_profile = Get_GeneBody_Profile(opt.bed_file, opt.gtf_file, opt.known_gene_list, opt.genic_partition, opt.fragment_size)
		
		profile_Up_genebody_Down_site(Upprofile, genebody_profile, Downprofile, opt.resolution, opt.genic_partition,
		opt.upstreamExtension, opt.downstreamExtension, opt.list_name, opt.condition, 'TSS_Body_TES')
	else:
		norm_Site_profile = Get_Site_Profile(opt.bed_file, opt.gtf_file, opt.known_gene_list, opt.type_site, opt.window_size,
		opt.fragment_size, opt.resolution, opt.upstreamExtension, opt.downstreamExtension)
	#### Second Site
		profile_plot_site(norm_Site_profile, opt.resolution, opt.upstreamExtension, opt.downstreamExtension, opt.list_name, opt.condition, opt.type_site)
		
		
		
if __name__ == "__main__":
	main(sys.argv)

