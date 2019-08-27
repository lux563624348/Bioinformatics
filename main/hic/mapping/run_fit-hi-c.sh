#!/bin/bash

EXEDIR=/home/zzeng/Software/fit-hi-c/bin

UPPER_BOUND=2000000
LOWER_BOUND=20000

for SAMPLE in stim_WT_CD8
do

	MARGINAL_CONTACTCOUNT_FILE=${SAMPLE}_marginal_contactCounts.txt.gz
	CONTACTCOUNT_FILE=${SAMPLE}_contactCounts.txt.gz
	ICE_BIAS_FILE=${SAMPLE}_ICE_bias.txt.gz

	if ! [ -f $MARGINAL_CONTACTCOUNT_FILE ]
	then
		gzip ${MARGINAL_CONTACTCOUNT_FILE%.txt.gz}.txt
	fi

	if ! [ -f $CONTACTCOUNT_FILE ]
        then
                gzip ${CONTACTCOUNT_FILE%.txt.gz}.txt
        fi

        if ! [ -f $ICE_BIAS_FILE ]
        then
                gzip ${ICE_BIAS_FILE%.txt.gz}.txt
        fi

	echo "python $EXEDIR/fit-hi-c.py -f $MARGINAL_CONTACTCOUNT_FILE -i $CONTACTCOUNT_FILE -t $ICE_BIAS_FILE -o ./ -l $SAMPLE -L $LOWER_BOUND -U $UPPER_BOUND"
	python $EXEDIR/fit-hi-c.py -f $MARGINAL_CONTACTCOUNT_FILE -i $CONTACTCOUNT_FILE -t $ICE_BIAS_FILE -o ./ -l $SAMPLE -L $LOWER_BOUND -U $UPPER_BOUND
	echo ""
	echo ""

done

echo "done"
