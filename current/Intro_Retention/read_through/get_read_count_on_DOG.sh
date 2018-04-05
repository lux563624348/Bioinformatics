#!/bin/bash

EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/IRTools/exploration/read_through

GTFDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Annotation/gtf_files
GTFFILE=hg19_genes.gtf

INDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/tophat_analysis
SAMPLE_LIST=(Naive_CD4_Rest Naive_CD4_Active_40min Naive_CD4_Active_150min Naive_CD4_Active_15h TCM_CD4_Rest TCM_CD4_Active_40min TCM_CD4_Active_150min TCM_CD4_Active_15h TEM_CD4_Rest TEM_CD4_Active_40min TEM_CD4_Active_150min TEM_CD4_Active_15h)

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/IRTools/exploration/read_through

READ_TYPE=single
LIBRARY_TYPE=fr-unstranded

DOG_EXTENSION=10000

for i in {0..11}
do

        SAMPLE=${SAMPLE_LIST[$i]}
        INPUTDIR=$INDIR/$SAMPLE


        OUTPUTFILE=${SAMPLE}_DOG_10K_read_count.txt

        echo "python get_read_count_on_DOG.py -b $INPUTDIR/${SAMPLE}.bam -p $READ_TYPE -s $LIBRARY_TYPE -g $GTFDIR/$GTFFILE -d $DOG_EXTENSION -o $OUTPUTDIR/$OUTPUTFILE"
        python get_read_count_on_DOG.py -b $INPUTDIR/${SAMPLE}.bam -p $READ_TYPE -s $LIBRARY_TYPE -g $GTFDIR/$GTFFILE -d $DOG_EXTENSION -o $OUTPUTDIR/$OUTPUTFILE
        echo ""
        echo ""

done

echo "done"
