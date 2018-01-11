#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

PATHTO=/SICER1.1
SICER=$PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 10 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["input dir"] ["ChIP bed file"] ["control bed file"] ["output dir"] ["species"] ["annotation dir"] ["annotation file"] ["gap size as multiply of window size"] ["effective genome fraction"] ["FDR"]  
    echo ""
    exit 1
fi

######################################################################
### Set up global variables input                                  ###
######################################################################

# The path to the data files.
DATADIR=$1
SAMPLEDIR=$DATADIR
CONTROLDIR=$DATADIR

# Output directory
OUTPUTDIR=$4
SUMMARY_DIR=$OUTPUTDIR
ISLANDS_BED_DIR=$OUTPUTDIR

# Input sample bed file
SAMPLEBED=$2
SAMPLE=${SAMPLEBED%.*}

# Input control bed file
CONTROLBED=$3
CONTROL=${CONTROLBED%.*}

# Species, for allowed species see GenomeData.py
SPECIES=$5

#Annotation file used in profile plotting
UCSCDIR=${6}
UCSCFILE=${7}

#Gap size in terms of multiply of window size
GAP=${8}

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=$9
EFFECTIVEGENOME=${EFFECTIVEGENOME:=0.74}

#EVALUE is the number of islands expected in random background. The E value is used for identification of candidate islands that exhibit clustering.
EVALUE=1000

#False discovery rate controlling significance
FDR=${10}


echo "###############################################################"
echo "######                    SICER v1.1                     ######"
echo "###############################################################"

echo "Input library directory: $DATADIR"
echo "ChIP library: $SAMPLEBED"
echo "Control library: $CONTROLBED"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Effective genome size as a fraction of the reference genome of $SPECIES: $EFFECTIVEGENOME"
echo "Evalue for identification of candidate islands that exhibit clustering: $EVALUE"
echo "False discovery rate controlling significance: $FDR"


######################################################################
### Remove redundant reads                                         ###
######################################################################

PVALUE=1e-5

# This file stores the preprocessed raw bed file.
FILTEREDSAMPLEBED=$SAMPLE-removed.bed
FILTEREDCONTROLBED=$CONTROL-removed.bed

if test -s $OUTPUTDIR/$FILTEREDSAMPLEBED; then
	echo ""
	echo ""
	echo " $OUTPUTDIR/$FILTEREDSAMPLEBED already exists, skipping this step ..."
	continue;
else
	echo " "
	echo " "
	echo "Preprocess the raw $SAMPLE file to remove redundancy with threshold $PVALUE..."
	echo "python $SICER/src/remove_redundant_reads_new.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $PVALUE -o $OUTPUTDIR/$FILTEREDSAMPLEBED"
	echo " "
	python $SICER/src/remove_redundant_reads_new.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -p $PVALUE -o $OUTPUTDIR/$FILTEREDSAMPLEBED
fi

if test -s $OUTPUTDIR/$FILTEREDCONTROLBED; then
	echo ""
	echo ""
	echo " $OUTPUTDIR/$FILTEREDCONTROLBED already exists, skipping this step ..."
	continue;
else
	echo " "
	echo " "
	echo "Preprocess the raw $CONTROL file to remove redundancy with threshold $PVALUE..."
	echo "python $SICER/src/remove_redundant_reads_new.py -s $SPECIES -b $CONTROLDIR/$CONTROLBED -t $PVALUE -o $OUTPUTDIR/$FILTEREDCONTROLBED"
	echo " "
	python $SICER/src/remove_redundant_reads_new.py -s $SPECIES -b $CONTROLDIR/$CONTROLBED -p $PVALUE -o $OUTPUTDIR/$FILTEREDCONTROLBED
fi

######################################################################
### Fragment size estimation                                       ###
######################################################################

RESOLUTION=5
FRAGMENT_SIZE=1

# Calculating tag correlation in the longest chromosome in the given species

CHROM=`python $SICER/lib/xs_lookup_chrom_length.py -s $SPECIES | grep $SPECIES | awk '{print $2}'`
let CHROM_LENGTH=`python $SICER/lib/xs_lookup_chrom_length.py -s $SPECIES | grep $SPECIES | awk '{print $3}'`
echo ""
echo "The longest chromosome in $SPECIES is $CHROM with length $CHROM_LENGTH"
echo ""

# Output file name
DATAFILE=$SAMPLE-$CHROM-tag-correlation

# Grep tags for longest chromosome
echo "Separate $CHROM tags to calculate fragment size ..."
grep ${CHROM}[[:space:]] $OUTPUTDIR/$FILTEREDCONTROLBED > temp1.bed
echo ""
echo "Total tag count in $CHROM is"
wc -l temp1.bed
echo ""

# Separate watson and crick reads
echo "Separate watson and crick reads in $CHROM ..."
grep [[:space:]]+ temp1.bed > temp_plus.bed
grep [[:space:]]- temp1.bed > temp_minus.bed
echo ""


# Make summary graph files
echo "Make summary graph files"
echo "Positive"
echo "python $SICER/lib/make_graph_file.py -f temp_plus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_plus.graph"
python $SICER/lib/make_graph_file.py -f temp_plus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_plus.graph
echo ""

echo ""
echo "Negative"
echo "python $SICER/lib/make_graph_file.py -f temp_minus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_minus.graph"
python $SICER/lib/make_graph_file.py -f temp_minus.bed -c $CHROM -l $CHROM_LENGTH -w $RESOLUTION -i 0 -o temp_minus.graph
echo ""

# Calculate cross correlation between watson summary graph and crick summary graph
echo ""
echo "Calculating correlation between watson and crick summary graphs for fragment size estimation"
echo "python $SICER/utility/calculate_cross_correlation_long_range.py -s $SPECIES -c $CHROM -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $OUTPUTDIR/$DATAFILE.txt"
let FRAGMENT_SIZE=`python $SICER/utility/calculate_cross_correlation_long_range.py -s $SPECIES -c $CHROM -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $OUTPUTDIR/$DATAFILE.txt | grep fragment | awk '{print $2}'`

if [ $FRAGMENT_SIZE -lt 100 ]; then
	echo ""
	echo "WARNING: CALCULATED FRAGMENT SIZE LESS THAN 100, USE 100 INSTEAD"
	echo ""
	FRAGMENT_SIZE=100
fi

if [ $FRAGMENT_SIZE -gt 300 ]; then
	echo ""
	echo "WARNING: CALCULATED FRAGMENT SIZE LARGER THAN 300, USE 300 INSTEAD"
	echo ""
	FRAGMENT_SIZE=300
fi

echo ""
echo "###############################################################"
echo "The estimated fragment size is: $FRAGMENT_SIZE"
echo "###############################################################"
echo ""

rm temp1.bed
rm temp_minus.bed
rm temp_minus.graph
rm temp_plus.bed
rm temp_plus.graph



######################################################################
### Window size calculation                                        ###
######################################################################

WINDOW_SIZE=1

echo ""
echo "Automatically calculate window size ..."
echo "python $SICER/extra/calculate-window-size.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -f $FRAGMENT_SIZE"
let WINDOW_SIZE=`python $SICER/extra/calculate-window-size.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -f $FRAGMENT_SIZE | grep Window | awk '{print $2}'`

let GAP_SIZE=$WINDOW_SIZE*$GAP

if [ $WINDOW_SIZE -lt 50 ]; then
	echo ""
	echo "WARNING: CALCULATED WINDOW SIZE LESS THAN 50"
	echo ""
fi

if [ $WINDOW_SIZE -gt 10000 ]; then
	echo ""
	echo "WARNING: CALCULATED WINDOW SIZE LARGER THAN 10000"
	echo ""
fi

echo ""
echo "################################################################"
echo "The selected Window Size is: $WINDOW_SIZE and Gap size $GAP_SIZE"
echo "################################################################"
echo ""




######################################################################
### SICER main module                                              ###
######################################################################


# This file stores the summary graph.  
SAMPLESUMMARY=$SAMPLE-W$WINDOW_SIZE.graph
CONTROLSUMMARY=$CONTROL-W$WINDOW_SIZE.graph
NORMALIZEDSAMPLESUMMARY=$SAMPLE-W$WINDOW_SIZE-normalized.graph
NORMALIZEDCONTROLSUMMARY=$CONTROL-W$WINDOW_SIZE-normalized.graph

# This file stores the histogram of window read-count.  
SAMPLESUMMARYHIST=$SAMPLE-W$WINDOW_SIZE.graphhist
CONTROLSUMMARYHIST=$CONTROL-W$WINDOW_SIZE.graphhist

# This file stores the sample summary graph in wig vstep format
NORMALIZEDSAMPLESUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-normalized.wig
# This file stores the control summary graph in wig vstep format 
NORMALIZEDCONTROLSUMMARYWIG=$CONTROL-W$WINDOW_SIZE-normalized.wig


#This file stores the candidate islands.
ISLAND=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE.scoreisland

#This file stores the summary of candidate islands, including chrom start end read-count_sample read-count-control pvalue, fold change and qvalue
ISLANDSIGNIFICANCEFILE=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-islands-summary

#This file stores the summary of significant islands identified with FDR criterion.
SIGNIFICANTISLANDS=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-islands-summary-FDR$FDR


#This file stores the significant islands in  "chr     start   end" format
ISLANDBED=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-island.bed
# This file stores the histogram of island scores
ISLANDSCOREHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR.islandscorehist
# This file stores the histogram of island lengths
ISLANDLENGTHHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR.islandlengthhist


# This file stores the island-filtered non-redundant raw reads 
ISLANDFILTEREDRAWBED=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.bed
# This file stores summary graph made by the island-filtered non-redundant raw reads 
ISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.graph
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads 
NORMALIZEDISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered-normalized.graph 
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads in wig vstep format
NORMALIZEDISLANDFILTEREDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered-normalized.wig


echo " "
echo " "
echo "Partion the genome in windows ..."
echo "Generate sample summary file ..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SAMPLESUMMARY 

if test -s $SUMMARY_DIR/$CONTROLSUMMARY; then
	echo ""
	echo ""
	echo " $SUMMARY_DIR/$CONTROLSUMMARY already exists, skipping this step ..."
	continue;
else
	echo ""
	echo "Generate control summary file ..."
	echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDCONTROLBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARY"
	python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDCONTROLBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$CONTROLSUMMARY 
fi


echo ""
echo ""
echo "Normalize sample summary graph by total non-redundant reads per million for $SAMPLE ..."
echo "python $SICER/src/normalize.py -i $SUMMARY_DIR/$SAMPLESUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARY"
python $SICER/src/normalize.py -i $SUMMARY_DIR/$SAMPLESUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARY

if test -s $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY; then
	echo ""
	echo ""
	echo " $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY already exists, skipping this step ..."
	continue;
else
	echo ""
	echo ""
	echo "Normalize control summary graph by total non-redundant reads per million for $CONTROL ..."
	echo "python $SICER/src/normalize.py -i $SUMMARY_DIR/$CONTROLSUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY"
	python $SICER/src/normalize.py -i $SUMMARY_DIR/$CONTROLSUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY
fi


echo ""
echo ""
echo "Convert the normalized sample summary graph into wig vstep format..."
echo "sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARY $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARY $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARYWIG $SAMPLE $WINDOW_SIZE


if test -f $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARYWIG; then
	echo ""
	echo ""
	echo " $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARYWIG already exists, skipping this step ..."
	continue;
else
	echo ""
	echo ""
	echo "Convert the normalized control summary graph into wig vstep format..."
	echo "sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARYWIG $CONTROL $WINDOW_SIZE"
	sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARYWIG $CONTROL $WINDOW_SIZE
fi

rm $SUMMARY_DIR/$NORMALIZEDSAMPLESUMMARY
rm $SUMMARY_DIR/$NORMALIZEDCONTROLSUMMARY


echo " "
echo " "
echo "Find candidate islands exhibiting clustering ..."
echo "python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE -f $ISLANDS_BED_DIR/$ISLAND"
python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SAMPLESUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE  -f $ISLANDS_BED_DIR/$ISLAND



echo ""
echo ""
echo "Calculate significance of candidate islands using the control library ..."
echo "python $SICER/src/associate_tags_with_chip_and_control_w_fc_q.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$ISLAND -f $FRAGMENT_SIZE -t $EFFECTIVEGENOME -o $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE"
python $SICER/src/associate_tags_with_chip_and_control_w_fc_q.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$ISLAND -f $FRAGMENT_SIZE -t $EFFECTIVEGENOME -o $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE



echo ""
echo ""
echo "Identify significant islands using FDR criterion ..."
echo "python $SICER/src/filter_islands_by_significance.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -p $FDR -c 7 -o $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS"
python $SICER/src/filter_islands_by_significance.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -p $FDR -c 7 -o $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS



echo ""
echo ""
echo "Convert island summary to island bed file of format chr start end ChIP-read-count"
echo "python $SICER/utility/convert_summary_to_bed.py -i $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS  -o  $ISLANDS_BED_DIR/$ISLANDBED"
python $SICER/utility/convert_summary_to_bed.py -i $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS  -o  $ISLANDS_BED_DIR/$ISLANDBED



echo ""
echo ""
echo "Filter reads with identified significant islands..."
echo "python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLANDBED  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED"
python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLANDBED  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED



echo ""
echo ""
echo "Make summary graph with filtered reads..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY



echo ""
echo ""
echo "Normalize summary graph with filtered reads for $SAMPLE by total island filtered reads per million..."
echo "python $SICER/src/normalize.py -i $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY -a 3 -t 1000000 -o $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY"
python $SICER/src/normalize.py -i $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY -a 3 -t 1000000 -o $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY



echo ""
echo ""
echo "Convert the summary graph made with the filtered reads into wig vstep format and normalize by total island-filtered read count per million..."
echo "sh $SICER/src/variableStep.sh $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARYWIG $SAMPLE-islandfiltered $WINDOW_SIZE

rm $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY
rm $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY



######################################################################
### Generate profiles around TSS, TES and genebody                 ###
######################################################################

EXEDIR=$SICER/extra
EXE=GenerateProfileAroundLocations.py

TYPE=TSS
NORMALIZATION=1
UPSTREAM=5000
DOWNSTREAM=5000
RESOLUTION=5
WINDOWSIZE=100
OUTFILE=${SAMPLE}_on_${TYPE}_profile_R${RESOLUTION}_W$WINDOWSIZE.txt

echo ""
echo ""
echo "Generating Profile around $TYPE..."
echo "python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE"

python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE

TYPE=TES
OUTFILE=${SAMPLE}_on_${TYPE}_profile_R${RESOLUTION}_W$WINDOWSIZE.txt

echo ""
echo ""
echo "Generating Profile around $TYPE..."
echo "python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE"

python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -n $NORMALIZATION -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE


EXE=GenerateProfileAroundRegions.py

TYPE=GENE
RESOLUTION=500
WINDOWSIZE=500
GENICPARTITION=40
let PLUSSHIFT=$FRAGMENT_SIZE/2
let MINUSSHIFT=$FRAGMENT_SIZE/2
OUTFILE=${SAMPLE}_on_${TYPE}_profile.txt

echo ""
echo ""
echo "Generating Profile around $TYPE..."
echo "python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -g $GENICPARTITION -p $PLUSSHIFT -m $MINUSSHIFT"

python $EXEDIR/$EXE  -k $UCSCDIR/$UCSCFILE -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -c $TYPE -o $OUTPUTDIR/$OUTFILE -s $SPECIES -u $UPSTREAM -d $DOWNSTREAM -r $RESOLUTION -w $WINDOWSIZE -g $GENICPARTITION -p $PLUSSHIFT -m $MINUSSHIFT



