#!/bin/bash

EXEDIR=/home/zzeng/Software
SPECIES=mm9

for SAMPLE in stim_WT_CD8
do

	echo "java -Xmx16g -jar $EXEDIR/juicer_tools_0.7.0.jar pre -n ${SAMPLE}_Juicebox_input.txt.gz ${SAMPLE}_Juicebox.hic $SPECIES"
	java -Xmx16g -jar $EXEDIR/juicer_tools_0.7.0.jar pre -n ${SAMPLE}_Juicebox_input.txt.gz ${SAMPLE}_Juicebox.hic $SPECIES
	echo ""
	echo ""

done

echo "done"
