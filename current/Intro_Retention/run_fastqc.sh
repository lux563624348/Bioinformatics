#!/bin/bash

INPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/raw_data
cd $INPUTDIR

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/fastqc

echo "fastqc -o $OUTPUTDIR $(ls *.gz | xargs)"
fastqc -o $OUTPUTDIR $(ls *.gz | xargs)
echo ""
echo ""

echo "done"
