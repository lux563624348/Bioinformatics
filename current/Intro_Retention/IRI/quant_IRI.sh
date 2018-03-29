#!/bin/bash

INDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/tophat_analysis
SAMPLE_LIST=(Naive_CD4_Rest Naive_CD4_Active_40min Naive_CD4_Active_150min Naive_CD4_Active_15h TCM_CD4_Rest TCM_CD4_Active_40min TCM_CD4_Active_150min TCM_CD4_Active_15h TEM_CD4_Rest TEM_CD4_Active_40min TEM_CD4_Active_150min TEM_CD4_Active_15h)

MAPDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/mappability
MAPFILE=hg19_wgEncodeCrgMapabilityAlign50mer.bigWig

for i in {0..11}
do

	SAMPLE=${SAMPLE_LIST[$i]}
	INPUTDIR=$INDIR/$SAMPLE

	echo "IRTools quant -q IRI -i $INPUTDIR/${SAMPLE}.bam -p single -s fr-unstranded -e hg19 -u $MAPDIR/$MAPFILE -n $SAMPLE"
	IRTools quant -q IRI -i $INPUTDIR/${SAMPLE}.bam -p single -s fr-unstranded -e hg19 -u $MAPDIR/$MAPFILE -n $SAMPLE
	echo ""
	echo ""

done

echo "done"
