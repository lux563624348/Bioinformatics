#!/bin/bash

EXEDIR=/home/data/SICER1.1/SICER/extra/tools/RNASeq
ALUDIR=/home/data/hg19/Annotation/RepElements
ALUSTRANSCRIPT=simplegenes_hg19_Alu_distribution_in_merged_transcript.pkl
ALUSINTRON=simplegenes_hg19_Alu_distribution_in_shared_intron.pkl
GENESDIR=/home/data/hg19/Annotation/NewVersion
GENES=hg19_EntrezID_filtered_samestrand_collisonremoved_simplegene.pkl
SPECIES=hg19
RNASEQDIR=/home/data/hg19/JunZhu/processed/Round2/RNASeq/newRUNOnRNASeq
RNASEQINTRON=Rest_R1_unique-on-EligibleEntrezGenes.dat_shared_intron_RPKMS.pkl
RNASEQEXON=Rest_R1_unique-on-EligibleEntrezGenes.dat_shared_exon_RPKMS.pkl
ENTREZIDS=IRI_high_and_down_entrez_ids.dat



echo "python $EXEDIR/Combine_Alu_RNASeq_intron_retention.py  -u $GENESDIR/$GENES -a $RNASEQDIR/$RNASEQINTRON -b $RNASEQDIR/$RNASEQEXON -d $ALUDIR/$ALUSINTRON -e $ALUDIR/$ALUSTRANSCRIPT -s hg19 -f $ENTREZIDS"

python $EXEDIR/Combine_Alu_RNASeq_intron_retention.py -u $GENESDIR/$GENES  -a $RNASEQDIR/$RNASEQINTRON -b $RNASEQDIR/$RNASEQEXON -d $ALUDIR/$ALUSINTRON -e $ALUDIR/$ALUSTRANSCRIPT -s hg19 -f $ENTREZIDS
