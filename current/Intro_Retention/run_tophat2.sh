#!/bin/bash

SICER=/home/lxiang/Software/SICER1.1/SICER
EXEDIR=$SICER/extra/tools/wiggle

THREADS=4

INPUTDIR=/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/raw_data
OUTPUTDIR=/home/lxiang/cloud_research/PengGroup/XLi/Raw_Data/Intro_Retention/Human/CD4

mkdir -p ${OUTPUTDIR}

SAMPLENAME_SET=(Naive_CD4_Rest Naive_CD4_Active_40min Naive_CD4_Active_150min Naive_CD4_Active_15h TCM_CD4_Rest TCM_CD4_Active_40min TCM_CD4_Active_150min TCM_CD4_Active_15h TEM_CD4_Rest TEM_CD4_Active_40min TEM_CD4_Active_150min TEM_CD4_Active_15h)

SPECIES=hg19
WINDOW_SIZE=10
FRAGMENT_SIZE=0

for i in {0..11}
do

	FASTQFILE=SRR48383$((21+$i)).fastq.gz
	cp ${INPUTDIR}/${FASTQFILE} ${OUTPUTDIR}/${SAMPLENAME_SET[$i]}.fastq.gz

	
done

echo "done"
