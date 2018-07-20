DIR=/home/wyang/Modules_WJ
bedDIR=/home/data/hg18/CD16/raw/set1

file1=GA2068-hg18-CD16-A1-RNAseq.bed
file2=GA2069-hg18-CD16-B1-RNAseq.bed

python $DIR/get_read_count_on_exons2_test.py -k knownGene_test.txt -b $bedDIR/$file1 -s hg18 -o Tags_on_$file1

#python $DIR/get_read_count_on_exons2.py -k /home/wyang/data/hg18/annotation/knownGene.txt -b $bedDIR/$file2 -s hg18 -o Tags_on_$file2

#sort -g -k 1 Tags_on_$file1 > Tags_on_$file1.sorted
#sort -g -k 1 Tags_on_$file2 > Tags_on_$file2.sorted

#cat Tags_on_$file1.sorted | awk '{print $2"\t"$1}' | uniq -u -f 1 | awk '{print $2"\t"$1}' >& Tags_on_$file1.unique

#cat Tags_on_$file2.sorted | awk '{print $2"\t"$1}' | uniq -u -f 1 | awk '{print $2"\t"$1}' >& Tags_on_$file2.unique


