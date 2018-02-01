#!/bin/bash

EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/paper_analysis/find_peak_summits/peaks_distribution

GTFDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/gtf_files
GTFFILE=mm9_genes.gtf

PROMOTER_UPSTREAM_EXTENSION=5000
TSS_REGION_LENGTH=2000

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/paper_analysis/find_peak_summits/peaks_distribution

INPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/paper_analysis/find_peak_summits/peaks_distribution

for PEAK_TYPE in WT_Tfh_Tcf1+_Ezh2+_peaks
do

	INPUTFILE=${PEAK_TYPE}.bed
	OUTPUTFILE=${PEAK_TYPE}_distribution.txt

	echo "python $EXEDIR/peaks_count_on_genic_region.py -i $INPUTDIR/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o $OUTPUTDIR/$OUTPUTFILE"
	python $EXEDIR/peaks_count_on_genic_region.py -i $INPUTDIR/$INPUTFILE -g $GTFDIR/$GTFFILE -u $PROMOTER_UPSTREAM_EXTENSION -t $TSS_REGION_LENGTH -o $OUTPUTDIR/$OUTPUTFILE
	echo ""
	echo ""

done

echo "done"
