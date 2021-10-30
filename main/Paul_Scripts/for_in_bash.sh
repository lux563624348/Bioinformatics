cd /slst/home/fanyh/Data/JGA/HISAT2_Results

xx=$(find -name *.sam)

for x in ${xx[*]}
do 
echo ${x}
samtools sort -o ${x:0:${#x}-4}_csorted.bam ${x} & # this is for old bash
done

