#!/bin/sh

DIR=/home/data/hg19/Annotation
#python Generate_Entrez_ucsc.py -u refFlat_hg19_EntrezID.ucsc -s hg19 -o refFlat_hg19_EntrezID.pkl
python ../Entrez.py -u $DIR/refFlat_hg19_EntrezID_filtered.ucsc -s hg19 -o hg19_EntrezID_filtered