#!/bin/bash
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


# This one should retire!



##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=/home/data/SICER1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 7 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["library file 1"] ["library file 2"] ["window size (bp)"] ["gap size (bp)"] ["E-value"] ["effective genome fraction for library 1"]  ["effective genome fraction for library 2"]
    echo ""
    exit 1
fi

TEMP=`expr $4 % $3`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
echo "#############################################"
echo "######           SICER v1.1            ######"
echo "#############################################"
# Input sample bed file
SAMPLEBED=$1
SAMPLE=${SAMPLEBED%.*}
echo "Library 1: $SAMPLEBED"

# Input control bed file
CONTROLBED=$2
CONTROL=${CONTROLBED%.*}
echo "Library 2: $CONTROLBED"

# Species, for allowed species see GenomeData.py
SPECIES=hg18
echo "Species: $SPECIES"

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME1=$6
echo "Effective genome size as a fraction of reference genome for $SAMPLEBED: $EFFECTIVEGENOME1"

EFFECTIVEGENOME2=$7
echo "Effective genome size as a fraction of reference genome for $CONTROLBED: $EFFECTIVEGENOME2"


# CHIPTHRESHOLD is the threshold is for redundancy allowed for reads. CHIPTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
CHIPTHRESHOLD=1
echo "Threshold for redundancy allowed for treated reads: $CHIPTHRESHOLD"

# CONTROLTHRESHOLD is the threshold is for redundancy allowed for reads. CONTROLTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
CONTROLTHRESHOLD=1
echo "Threshold for redundancy allowed for control reads: $CONTROLTHRESHOLD"

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$3
echo "Window size: $WINDOW_SIZE bps"

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=150
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"

#GAP_SIZE is in base pairs.
GAP_SIZE=$4
echo "Gap size: $GAP_SIZE bps"

#EVALUE is the number of islands expected in random background. The E value is used for identification of candidate islands that exhibit clustering.
EVALUE=$5
echo "Evalue for identification of candidate islands that exhibit clustering: $EVALUE"


#================================================================================
#############################################
# ######  SET UP DIRECTORY STRUCTURE ###### #
#############################################
# If data files are not in the current directory, replace this with
# the path to the data files.
DATADIR=.
SAMPLEDIR=/home/data/hg18/CD8/2010data/092010
CONTROLDIR=/home/data/hg18/CD8/2008data
# If You want the output files not in the current directory, replace this.
OUTPUTDIR=.
SUMMARY_DIR=$OUTPUTDIR
ISLANDS_BED_DIR=$OUTPUTDIR

#================================================================================
#############################################
# ###### DEFINITION OF OUTPUT FILES  ###### #
#############################################

# This file stores the preprocessed raw bed file.
FILTEREDSAMPLEBED=$SAMPLE-${CHIPTHRESHOLD}-removed.bed
FILTEREDCONTROLBED=$CONTROL-${CONTROLTHRESHOLD}-removed.bed

# This file stores the summary graph.  
SAMPLESUMMARY=$SAMPLE-W$WINDOW_SIZE.graph
CONTROLSUMMARY=$CONTROL-W$WINDOW_SIZE.graph

# This file stores the histogram of window read-count.  
SAMPLESUMMARYHIST=$SAMPLE-W$WINDOW_SIZE.graphhist
CONTROLSUMMARYHIST=$CONTROL-W$WINDOW_SIZE.graphhist

#This file stores the candidate islands.
TREATEDISLAND=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE.island
CONTROLISLAND=$CONTROL-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE.island

# This file stores the histogram of island scores
SAMPLEISLANDSCOREHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE.islandscorehist
# This file stores the histogram of island lengths
SAMPLEISLANDLENGTHHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE.islandlengthhist
# This file stores the histogram of island scores
CONTROLISLANDSCOREHIST=$CONTROL-W$WINDOW_SIZE-G$GAP_SIZE.islandscorehist
# This file stores the histogram of island lengths
CONTROLISLANDLENGTHHIST=$CONTROL-W$WINDOW_SIZE-G$GAP_SIZE.islandlengthhist


UNIONISLAND=$SAMPLE-and-$CONTROL-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-union.island


#This file stores the summary of candidate islands, including tags counts, 
# pvalue, fold change and BH-corrected p-value
COMPARISONSUMMARY=$SAMPLE-vs-$CONTROL-W$WINDOW_SIZE-G$GAP_SIZE-comparison-summary


echo " "
echo " "
echo "Preprocess the raw $SAMPLE file to remove redundancy with threshold $CHIPTHRESHOLD..."
echo "python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $CHIPTHRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED"
python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $CHIPTHRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED

echo " "
echo " "
echo "Preprocess the raw $CONTROL file to remove redundancy with threshold $CONTROLTHRESHOLD..."
echo "python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $CONTROLDIR/$CONTROL.bed -t $CONTROLTHRESHOLD -o $OUTPUTDIR/$FILTEREDCONTROLBED"
python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $CONTROLDIR/$CONTROL.bed -t $CONTROLTHRESHOLD -o $OUTPUTDIR/$FILTEREDCONTROLBED


echo " "
echo " "
echo "Partion the genome in windows ..."
echo "Generate summary files for $SAMPLE..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARY 





# Generating window read-count statistics
echo " "
echo " "
echo "Generate window read count histograms for $SAMPLE ... "
echo "python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARYHIST"
python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARYHIST

echo " "
echo " "
echo "Generate summary files for $CONTROL..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDCONTROLBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDCONTROLBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARY 


# Generating window read-count statistics
echo " "
echo " "
echo "Generate window read count histograms for $CONTROL ... "
echo "python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARYHIST"
python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARYHIST




echo " "
echo " "
echo "Find significant islands with E-value $EVALUE for $SAMPLE..."
echo "python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME1 -e $EVALUE -f $ISLANDS_BED_DIR/$TREATEDISLAND"
python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME1 -e $EVALUE  -f $ISLANDS_BED_DIR/$TREATEDISLAND

echo " "
echo " "
echo "Find significant islands with E-value $EVALUE for $CONTROL..."
echo "python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME2 -e $EVALUE  -f $ISLANDS_BED_DIR/$CONTROLISLAND"
python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME2 -e $EVALUE  -f $ISLANDS_BED_DIR/$CONTROLISLAND


#Generating island statistics, optional
echo " "
echo " "
echo "Get island statistics for $SAMPLE..."
echo "python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -i $ISLANDS_BED_DIR/$TREATEDISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$SAMPLEISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$SAMPLEISLANDLENGTHHIST"
python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -i $ISLANDS_BED_DIR/$TREATEDISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$SAMPLEISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$SAMPLEISLANDLENGTHHIST


#Generating island statistics, optional
echo " "
echo " "
echo "Get island statistics for $CONTROL..."
echo "python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -i $ISLANDS_BED_DIR/$CONTROLISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$CONTROLISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$CONTROLISLANDLENGTHHIST"
python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$CONTROLSUMMARY -i $ISLANDS_BED_DIR/$CONTROLISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$CONTROLISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$CONTROLISLANDLENGTHHIST



echo ""
echo ""
echo "Merge the two identified sets of significant islands..."
echo "python $SICER/src/find_union_islands.py -s $SPECIES -a $ISLANDS_BED_DIR/$CONTROLISLAND -b $ISLANDS_BED_DIR/$TREATEDISLAND -o $ISLANDS_BED_DIR/$UNIONISLAND"
python $SICER/src/find_union_islands.py -s $SPECIES -a $ISLANDS_BED_DIR/$CONTROLISLAND -b $ISLANDS_BED_DIR/$TREATEDISLAND -o $ISLANDS_BED_DIR/$UNIONISLAND

echo ""
echo ""
echo "Find read counts from the two libraries on the merged islands and calculate correlation"
echo "python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$COMPARISONSUMMARY"

python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$COMPARISONSUMMARY

echo ""
echo ""
echo "Done!"
