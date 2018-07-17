#!/bin/bash

EXEDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/Ezh2_ChIPSeq_Sep2016/SICER_1E-8/processed/profiles_around_Tfh_Ezh2_Tcf1_common_peaks

INDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/Ezh2_ChIPSeq_Sep2016/SICER_1E-8

OUTPUTDIR=/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/Ezh2_ChIPSeq_Sep2016/SICER_1E-8/processed/profiles_around_Tfh_Ezh2_Tcf1_common_peaks

SITETYPE_LIST=(new_Ezh2_Tcf1_common_sites Ezh2_Tcf1_common_sites)
SITEDIR_LIST=(/home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/Ezh2_ChIPSeq_Sep2016/SICER_1E-8/processed/overlap_with_Tcf1_peaks /home/zzeng/cloud_research/PengGroup/ZZeng/Data/Haihui/Ezh2/ChIP-seq/SICER_1E-8/processed/overlap_with_Tcf1_peaks)
SITEFILE_LIST=(new_Tfh_Ezh2_Tcf1_common_peaks.bed Tfh_Ezh2_Tcf1_common_peaks.bed)

NORMALIZATION=1
FRAGMENTSIZE=150
UPSTREAMEXTENSION=5000
DOWNSTREAMEXTENSION=5000
WINDOWSIZE=250
RESOLUTION=50
GENIC_PARTION=20
DRAW_PROFILE_PLOT=no

for i in {0..1}
do

SITETYPE=${SITETYPE_LIST[$i]}
SITEDIR=${SITEDIR_LIST[$i]}
SITEFILE=${SITEFILE_LIST[$i]}

GAP=600

for CELLTYPE in WT_Tfh DKO_Tfh
do
    for HISMODS in Ezh2
    do

	READDIR=$INDIR/${CELLTYPE}_${HISMODS}
	READFILE=${CELLTYPE}_${HISMODS}-W200-G${GAP}-FDR0.00000001-islandfiltered.bed
	OUTFILE=${CELLTYPE}_${HISMODS}_on_${SITETYPE}.txt

	echo "python $EXEDIR/generate_profile_around_sites.py -b $READDIR/$READFILE -s $SITEDIR/$SITEFILE -n $NORMALIZATION -f $FRAGMENTSIZE -r $RESOLUTION -u $UPSTREAMEXTENSION -d $DOWNSTREAMEXTENSION -w $WINDOWSIZE -p $GENIC_PARTION -y $DRAW_PROFILE_PLOT -o $OUTPUTDIR/$OUTFILE"
	python $EXEDIR/generate_profile_around_sites.py -b $READDIR/$READFILE -s $SITEDIR/$SITEFILE -n $NORMALIZATION -f $FRAGMENTSIZE -r $RESOLUTION -u $UPSTREAMEXTENSION -d $DOWNSTREAMEXTENSION -w $WINDOWSIZE -p $GENIC_PARTION -y $DRAW_PROFILE_PLOT -o $OUTPUTDIR/$OUTFILE
	echo ""
	echo ""
    
    done
done

done


echo "done"
