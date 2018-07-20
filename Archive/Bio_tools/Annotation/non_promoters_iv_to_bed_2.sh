#!/bin/bash
EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/gene_iv/mm9

GTFDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/gtf_files
GTFFILE=mm9_genes.gtf

UPSTREAM_EXTENSION=1000
DOWNSTREAM_EXTENSION=0

REGION=NonPromoter

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/gene_iv/mm9
OUTPUTFILE=gene_non-promoter_1k_iv.bed

echo "python $EXEDIR/promoters_iv_to_bed.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE"
python $EXEDIR/promoters_iv_to_bed.py -g $GTFDIR/$GTFFILE -u $UPSTREAM_EXTENSION -d $DOWNSTREAM_EXTENSION -r $REGION -o $OUTPUTDIR/$OUTPUTFILE
echo ""
echo ""

echo "done"
