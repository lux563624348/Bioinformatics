#!/bin/bash
# Copyright (c) 2011 The George Washington University
# Authors: Weiqun Peng
#

 
sh ./SICER-df.sh GA1564_hg18_CD16-H3K4me1-A-n500000.bed GA1566_hg18_CD16-IgG-A-n500000.bed GA1568_hg18_CD16-H3K4me1-B-n500000.bed GA1570_hg18_CD16-IgG-B-n500000.bed 200 400 .001 .001
echo "sh ./SICER-df.sh GA1564_hg18_CD16-H3K4me1-A-n500000.bed GA1566_hg18_CD16-IgG-A-n500000.bed GA1568_hg18_CD16-H3K4me1-B-n500000.bed GA1570_hg18_CD16-IgG-B-n500000.bed 200 400 .001 .001"

#["KO bed file"] ["KO control file"] ["WT bed file"]  ["WT control file"] ["window size (bp)"] ["gap size (bp)"] ["FDR for KO vs KOCONTROL or WT vs WTCONTROL"] ["FDR for WT vs KO"]