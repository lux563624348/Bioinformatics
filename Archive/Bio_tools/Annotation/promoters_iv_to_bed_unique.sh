#!/bin/bash
EXEDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm9

GTFDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gtf_files
GTFFILE=mm9_genes.gtf

UPSTREAM_EXTENSION=0
DOWNSTREAM_EXTENSION=1000

REGION=Promoter

OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Annotation/gene_iv/mm9
OUTPUTFILE=gene_genebody_iv_unique.bed

echo "python $EXEDIR/promoters_iv_to_bed_unique.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE"
python $EXEDIR/promoters_iv_to_bed_unique.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE
echo ""
echo ""

echo "done"
