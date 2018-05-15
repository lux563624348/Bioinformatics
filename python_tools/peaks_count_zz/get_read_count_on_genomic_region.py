import HTSeq
import sys
import collections
from optparse import OptionParser


def main(argv):
    parser = OptionParser()

    parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedFile", metavar="<file>",
                      help="ChIP seq read file")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size",
                      help="fragment_size determins the shift (half of fragment_size of ChIP-seq read position, in bps",
                      metavar="<int>")
    parser.add_option("-g", "--genes_GTF_file", action="store", type="string", dest="genesfile", metavar="<file>",
                      help="GTF genes annotation file")
    parser.add_option("-k", "--gene_list", action="store", type="string", dest="gene_list", metavar="<file>",
                      help="gene list")
    parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'PromoterGenebody' or 'ExonicRegion'", action="store",
                      type="string", dest="region_type", metavar="<str>", help="region to count tags in")
    parser.add_option("-u", "--promoter_upstream_extension", action="store", type="int",
                      dest="promoter_upstream_extension", help="upstream extension of promoter region from TSS",
                      metavar="<int>")
    parser.add_option("-d", "--promoter_downstream_extension", action="store", type="int",
                      dest="promoter_downstream_extension", help="downstream extension of promoter region from TSS",
                      metavar="<int>")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>",
                      help="output file name for genes and tag numbers")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    if opt.gene_list:
        gene_list = []
        f = open(opt.gene_list, "r")
        for line in f:
            line = line.strip()
            sline = line.split()
            gene_list.append(sline[0])

    fragment_size = opt.fragment_size

    if opt.promoter_upstream_extension is not None:
        UpstreamExtension = opt.promoter_upstream_extension
    if opt.promoter_downstream_extension is not None:
        DownstreamExtension = opt.promoter_downstream_extension

    if opt.region_type == "Promoter":
        promoters = HTSeq.GenomicArrayOfSets("auto", stranded=False)

        counts = collections.Counter()
        promoter_region_length = collections.Counter()
        total_count = 0

        if opt.gene_list:
            gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
            for feature in gtffile:
                if feature.type == "transcript_region" and feature.attr["gene_id"] in gene_list:
                    transcript_region_iv = feature.iv
                    if transcript_region_iv.strand == "+":
                        if transcript_region_iv.start_d - UpstreamExtension < 0:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                transcript_region_iv.start_d + DownstreamExtension + 1,
                                                                transcript_region_iv.strand)
                        else:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                transcript_region_iv.start_d - UpstreamExtension,
                                                                transcript_region_iv.start_d + DownstreamExtension + 1,
                                                                transcript_region_iv.strand)
                    elif transcript_region_iv.strand == "-":
                        if transcript_region_iv.start_d - DownstreamExtension < 0:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                transcript_region_iv.start_d + UpstreamExtension + 1,
                                                                transcript_region_iv.strand)
                        else:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                transcript_region_iv.start_d - DownstreamExtension,
                                                                transcript_region_iv.start_d + UpstreamExtension + 1,
                                                                transcript_region_iv.strand)
                    gene_id = feature.attr["gene_id"]
                    promoters[promoter_iv] += gene_id
                    if gene_id not in counts.keys():
                        counts[gene_id] = 0
        else:
            gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
            for feature in gtffile:
                if feature.type == "transcript_region":
                    transcript_region_iv = feature.iv
                    if transcript_region_iv.strand == "+":
                        if transcript_region_iv.start_d - UpstreamExtension < 0:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                transcript_region_iv.start_d + DownstreamExtension + 1,
                                                                transcript_region_iv.strand)
                        else:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                transcript_region_iv.start_d - UpstreamExtension,
                                                                transcript_region_iv.start_d + DownstreamExtension + 1,
                                                                transcript_region_iv.strand)
                    elif transcript_region_iv.strand == "-":
                        if transcript_region_iv.start_d - DownstreamExtension < 0:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                transcript_region_iv.start_d + UpstreamExtension + 1,
                                                                transcript_region_iv.strand)
                        else:
                            promoter_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                transcript_region_iv.start_d - DownstreamExtension,
                                                                transcript_region_iv.start_d + UpstreamExtension + 1,
                                                                transcript_region_iv.strand)
                    gene_id = feature.attr["gene_id"]
                    promoters[promoter_iv] += gene_id
                    if gene_id not in counts.keys():
                        counts[gene_id] = 0

        for iv, step_set in promoters.steps():
            if len(step_set) != 0:
                for gene_id in list(step_set):
                    promoter_region_length[gene_id] += iv.length

        bedfile = HTSeq.BED_Reader(opt.bedFile)
        for alt in bedfile:
            total_count += 1
            if alt.iv.strand == "+":
                alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2, ".")
            elif alt.iv.strand == "-":
                alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2, ".")

            set_pos = promoters[alt_pos]
            if len(set_pos) == 1:
                gene_id = list(set_pos)[0]
                counts[gene_id] += 1

        f = open(opt.outfile, "w")
        outline = "#GeneName" + '\t' + "Read_Count" + '\t' + "RPKM" + '\t' + "RPKM_pseudo_count" + '\n'
        f.write(outline)

        for gene_id in sorted(counts.keys()):
            RPKM = counts[gene_id] / (promoter_region_length[gene_id] / 1000.0) / (total_count / 1000000.0)
            RPKM_pseudo_count = (counts[gene_id] + 1) / (promoter_region_length[gene_id] / 1000.0) / (
            total_count / 1000000.0)
            outline = gene_id + "\t" + str(counts[gene_id]) + "\t" + str(RPKM) + "\t" + str(RPKM_pseudo_count) + "\n"
            f.write(outline)

        f.close()

    elif opt.region_type == "GeneBody":
        pass
    elif opt.region_type == "PromoterGenebody":
        promoter_plus_genebody = HTSeq.GenomicArrayOfSets("auto", stranded=False)

        counts = collections.Counter()
        promoter_plus_genebody_region_length = collections.Counter()
        total_count = 0

        if opt.gene_list:
            gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
            for feature in gtffile:
                if feature.type == "transcript_region" and feature.attr["gene_id"] in gene_list:
                    transcript_region_iv = feature.iv
                    if transcript_region_iv.strand == "+":
                        if transcript_region_iv.start - UpstreamExtension < 0:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                              transcript_region_iv.end + DownstreamExtension,
                                                                              transcript_region_iv.strand)
                        else:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                              transcript_region_iv.start - UpstreamExtension,
                                                                              transcript_region_iv.end + DownstreamExtension,
                                                                              transcript_region_iv.strand)
                    elif transcript_region_iv.strand == "-":
                        if transcript_region_iv.start - DownstreamExtension < 0:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                              transcript_region_iv.end + UpstreamExtension,
                                                                              transcript_region_iv.strand)
                        else:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                              transcript_region_iv.start - DownstreamExtension,
                                                                              transcript_region_iv.end + UpstreamExtension,
                                                                              transcript_region_iv.strand)
                    gene_id = feature.attr["gene_id"]
                    promoter_plus_genebody[promoter_plus_genebody_iv] += gene_id
                    if gene_id not in counts.keys():
                        counts[gene_id] = 0
        else:
            gtffile = HTSeq.GFF_Reader(opt.genesfile, end_included=True)
            for feature in gtffile:
                if feature.type == "transcript_region":
                    transcript_region_iv = feature.iv
                    if transcript_region_iv.strand == "+":
                        if transcript_region_iv.start - UpstreamExtension < 0:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                              transcript_region_iv.end + DownstreamExtension,
                                                                              transcript_region_iv.strand)
                        else:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                              transcript_region_iv.start - UpstreamExtension,
                                                                              transcript_region_iv.end + DownstreamExtension,
                                                                              transcript_region_iv.strand)
                    elif transcript_region_iv.strand == "-":
                        if transcript_region_iv.start - DownstreamExtension < 0:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom, 0,
                                                                              transcript_region_iv.end + UpstreamExtension,
                                                                              transcript_region_iv.strand)
                        else:
                            promoter_plus_genebody_iv = HTSeq.GenomicInterval(transcript_region_iv.chrom,
                                                                              transcript_region_iv.start - DownstreamExtension,
                                                                              transcript_region_iv.end + UpstreamExtension,
                                                                              transcript_region_iv.strand)
                    gene_id = feature.attr["gene_id"]
                    promoter_plus_genebody[promoter_plus_genebody_iv] += gene_id
                    if gene_id not in counts.keys():
                        counts[gene_id] = 0

        for iv, step_set in promoter_plus_genebody.steps():
            if len(step_set) != 0:
                for gene_id in list(step_set):
                    promoter_plus_genebody_region_length[gene_id] += iv.length

        bedfile = HTSeq.BED_Reader(opt.bedFile)
        for alt in bedfile:
            total_count += 1
            if alt.iv.strand == "+":
                alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2, ".")
            elif alt.iv.strand == "-":
                alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2, ".")

            set_pos = promoter_plus_genebody[alt_pos]
            if len(set_pos) == 1:
                gene_id = list(set_pos)[0]
                counts[gene_id] += 1

        f = open(opt.outfile, "w")
        outline = "gene_id" + '\t' + "read_Count" + '\t' + "RPKM" + '\t' + "RPKM_with_pseudo_count" + '\n'
        f.write(outline)

        for gene_id in sorted(counts.keys()):
            RPKM = counts[gene_id] / (promoter_plus_genebody_region_length[gene_id] / 1000.0) / (
            total_count / 1000000.0)
            RPKM_pseudo_count = (counts[gene_id] + 1) / (promoter_plus_genebody_region_length[gene_id] / 1000.0) / (
            total_count / 1000000.0)
            outline = gene_id + "\t" + str(counts[gene_id]) + "\t" + str(RPKM) + "\t" + str(RPKM_pseudo_count) + "\n"
            f.write(outline)

        f.close()


if __name__ == "__main__":
    main(sys.argv)
