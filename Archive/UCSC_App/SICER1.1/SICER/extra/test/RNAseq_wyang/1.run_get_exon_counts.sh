exeDIR=/home/wyang/Modules_WJ
bedDIR=/home/data/hg18/CD16/raw/set1

file1=Set1_tags_on_Exons.txt
## get tags on exons of each gene
python $exeDIR/get_read_count_on_exons2.py -k /home/wyang/data/hg18/annotation/knownGene.txt -b $bedDIR/$file1 -s hg18 -o Tags_on_$file1

## keep unique genes
sort -g -k 1 $file1 > $file1.sorted
awk '!x[$1]++' $file1.sorted > $file1.unique 
