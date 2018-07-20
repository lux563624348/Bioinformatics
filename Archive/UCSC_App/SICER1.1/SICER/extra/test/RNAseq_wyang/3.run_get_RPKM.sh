DIR=/home/wyang/Modules_WJ

bedDIR1=/home/data/hg18/CD16/raw/set1
file1=GA2068-hg18-CD16-A1-RNAseq.bed

python $DIR/get_RPKM.py -k /home/wyang/data/hg18/annotation/knownGene.txt -b $bedDIR1/$file1 -s hg18 -o RPKM_on_$file1



