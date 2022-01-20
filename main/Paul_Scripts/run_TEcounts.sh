WorkDir=/public/slst/home/fanyh/Data/JGA/HISAT2_Results
TE_GTFFILE=/public/slst/home/fanyh/Annotation/GRCh38_GENCODE_rmsk_TE.gtf
GENE_GTFFILE=/public/slst/home/fanyh/Annotation/Homo_sapiens.GRCh38.105.gtf
cd ${WorkDir}

OutDir=/public/slst/home/fanyh/Data/JGA/TE_Results

bamFiles=$(find -name "*.bam")

start_date=`date`
i=0
for bam in ${bamFiles[*]};
do
    if [ ! -f ${OutDir}/${bam:2:-22}_count.cntTable ];then
        echo ${bam:2:-22}
        TEcount --sortByPos \
            --format BAM \
            --mode multi \
            -b ${bam} \
            --project ${OutDir}/${bam:2:-22}_count \
            --GTF ${GENE_GTFFILE} \
            --verbose 0 \
            --TE ${TE_GTFFILE} > ${OutDir}/log/${bam:2:-22}_count.log  &
        i=$(expr ${i} + 1)
    fi
    if ( expr ${i} % 10  == 9 )
    then
        #break
        sleep 1500s
    fi
done


end_date=`date`
echo "TEcount Start at ${start_date}" #| mail -s "Project Finished at ${end_date}" lux@gwu.edu

