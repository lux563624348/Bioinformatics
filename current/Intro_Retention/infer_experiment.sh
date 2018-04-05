#!/bin/bash

REFGENE_BED_DIR=/home/zzeng/Software/RSeQC-2.6.2/gene_model
REFGENE_BED_FILE=hg19_UCSC_knownGene.bed

BAM_DIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/intron_retention/Human/CD4/tophat_analysis/Naive_CD4_Active_15h
BAM_FILE=Naive_CD4_Active_15h.bam

echo "infer_experiment.py -r $REFGENE_BED_DIR/$REFGENE_BED_FILE -i $BAM_DIR/$BAM_FILE"
infer_experiment.py -r $REFGENE_BED_DIR/$REFGENE_BED_FILE -i $BAM_DIR/$BAM_FILE

echo "done"
