#!/bin/bash
EXEDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm10

GTFDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gtf_files
GTFFILE=mm10_genes.gtf


UPSTREAM_EXTENSION=2000
DOWNSTREAM_EXTENSION=2000

REGION=Promoter

OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm10
OUTPUTFILE=gene_promoter_2k_iv_unique_mm10.bed

echo "python $EXEDIR/promoters_iv_to_bed_unique.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE"
python $EXEDIR/promoters_iv_to_bed_unique.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE
echo ""
echo ""

echo "done"
