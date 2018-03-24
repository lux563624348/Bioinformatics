#!/bin/bash

EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/merge_peaks

BEDDIR_C1=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/WT_CD4_Ezh2
BEDFILE_C1=WT_CD4_Ezh2-1-removed.bed

BEDDIR_C2=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/WT_Tfh_Ezh2
BEDFILE_C2=WT_Tfh_Ezh2-1-removed.bed

PEAKSDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/merge_peaks
PEAKSFILE=WT_CD4_Tfh_merge_Ezh2.bed

FRAGMENTSIZE=150

SUMMARYDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/merge_peaks
SUMMARYFILE=read_count_on_WT_CD4_WT_Tfh_merge_Ezh2_peaks.summary

echo "python $EXEDIR/get_two_conditions_read_count_on_sites.py -a $BEDDIR_C1/$BEDFILE_C1 -b $BEDDIR_C2/$BEDFILE_C2 -p $PEAKSDIR/$PEAKSFILE -f $FRAGMENTSIZE -s $SUMMARYDIR/$SUMMARYFILE"
python $EXEDIR/get_two_conditions_read_count_on_sites.py -a $BEDDIR_C1/$BEDFILE_C1 -b $BEDDIR_C2/$BEDFILE_C2 -p $PEAKSDIR/$PEAKSFILE -f $FRAGMENTSIZE -s $SUMMARYDIR/$SUMMARYFILE
echo ""
echo ""


BEDDIR_C1=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/WT_CD8_Ezh2
BEDFILE_C1=WT_CD8_Ezh2-1-removed.bed

BEDDIR_C2=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/WT_CD8effector_Ezh2
BEDFILE_C2=WT_CD8effector_Ezh2-1-removed.bed

PEAKSDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/merge_peaks
PEAKSFILE=WT_CD8_CD8effector_merge_Ezh2.bed

FRAGMENTSIZE=100

SUMMARYDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER/processed/merge_peaks
SUMMARYFILE=read_count_on_WT_CD8_WT_CD8effector_merge_Ezh2_peaks.summary

echo "python $EXEDIR/get_two_conditions_read_count_on_sites.py -a $BEDDIR_C1/$BEDFILE_C1 -b $BEDDIR_C2/$BEDFILE_C2 -p $PEAKSDIR/$PEAKSFILE -f $FRAGMENTSIZE -s $SUMMARYDIR/$SUMMARYFILE"
python $EXEDIR/get_two_conditions_read_count_on_sites.py -a $BEDDIR_C1/$BEDFILE_C1 -b $BEDDIR_C2/$BEDFILE_C2 -p $PEAKSDIR/$PEAKSFILE -f $FRAGMENTSIZE -s $SUMMARYDIR/$SUMMARYFILE
echo ""
echo ""

echo "done"
