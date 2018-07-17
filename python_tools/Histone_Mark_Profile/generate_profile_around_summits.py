__author__ = 'Zhouhao Zeng'

import HTSeq
import sys
from optparse import OptionParser
import numpy


def getSummitProfile(ga, summit_pos_set, window_size, resolution, UpstreamExtension, DownstreamExtension):
    upstream_num_points = UpstreamExtension / resolution
    downstream_num_points = DownstreamExtension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1

    profile = numpy.zeros(total_num_points, numpy.int)

    for summit_pos in summit_pos_set:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            index_pos = summit_pos.pos + (index - upstream_num_points) * resolution
            summit_pos_window_iv = HTSeq.GenomicInterval(summit_pos.chrom, index_pos - window_size / 2,
                                                         index_pos + window_size / 2)

            for step_iv, step_count in ga[summit_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[index] += count_in_window
            index += 1

    return profile


def main(argv):
    desc = """This is a template for the analysis of aggretated tag distribution with respect to a set of points, such as the TSSs of known genes, with one profile from each strand."""
    parser = OptionParser(description=desc)
    parser.add_option("-s", "--summits_file", action="store", type="string",
                      dest="summitsfile", metavar="<file>", help="summits bed file")
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bed_file", help="file with tags in bed format", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="outfile name", metavar="<file>")
    parser.add_option("-n", "--normalization", action="store", type="float", dest="norm",
                      help="additional normalization in addition to number of sites, number of reads per million and window_size per 1K")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size",
                      help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps",
                      metavar="<int>")
    parser.add_option("-u", "--UpstreamExtension", action="store", type="int",
                      dest="upstreamExtension", help="UpstreamExtension", metavar="<int>")
    parser.add_option("-d", "--DownstreamExtension", action="store", type="int",
                      dest="downstreamExtension", help="DownstreamExtension", metavar="<int>")
    parser.add_option("-w", "--WindowSize", action="store", type="int",
                      dest="window_size",
                      help="window size for averaging. When window size > resolution, there is smoothing",
                      metavar="<int>")
    parser.add_option("-r", "--resolution", action="store", type="int",
                      dest="resolution", help="resolution of the upstream and downstream profile, eg, 5",
                      metavar="<int>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 18:
        parser.print_help()
        sys.exit(1)

    fragment_size = opt.fragment_size
    window_size = opt.window_size
    UpstreamExtension = opt.upstreamExtension
    DownstreamExtension = opt.downstreamExtension
    resolution = opt.resolution

    print "Upstream extension: %i" % UpstreamExtension
    print "Downstream extension: %i" % DownstreamExtension
    print "Scanning window size: %i" % window_size
    print "Scanning resolution: %i" % resolution

    num_tags = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    bedfile = HTSeq.BED_Reader(opt.bed_file)
    for alt in bedfile:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        ga[alt_pos] += 1
        num_tags += 1

    num_summits = 0
    summit_pos_set = set()
    summitsfile = HTSeq.BED_Reader(opt.summitsfile)
    for alt in summitsfile:
        summit_pos_set.add(alt.iv.start_d_as_pos)
        num_summits += 1

    profile = getSummitProfile(ga, summit_pos_set, window_size, resolution, UpstreamExtension, DownstreamExtension)

    normalization = num_tags / 1000000.0
    normalization *= num_summits
    normalization *= window_size / 1000.0
    normalization *= opt.norm

    print "Number of locations: %i" % num_summits
    print "Number of reads: %i" % num_tags
    print "Normalization is by total number of reads per million. normalization = %f" % normalization

    f = open(opt.outfile, "w")
    xValues = numpy.arange(-UpstreamExtension, DownstreamExtension + 1, resolution)
    normalized_profile = profile / normalization
    for index in range(len(xValues)):
        outline = str(xValues[index]) + "\t" + str(normalized_profile[index]) + "\n"
        f.write(outline)
    f.close()


if __name__ == "__main__":
    main(sys.argv)
