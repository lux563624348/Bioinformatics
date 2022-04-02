#WorkDir=/public/slst/home/fanyh/Data/JGA/Unique_Align/HISAT2_Results/
WorkDir=/slst/home/fanyh/Data/JGA/Multi_Align/HISAT2_Results
TE_GTFFILE=/public/slst/home/fanyh/Annotation/GRCh38_GENCODE_rmsk_TE.gtf
GENE_GTFFILE=/public/slst/home/fanyh/Annotation/Homo_sapiens.GRCh38.105.gtf

cd ${WorkDir}
bamFiles=$(find -name "*.bam")

start_date=`date`
i=0
LenAffix=".sorted.bam"


for bam in ${bamFiles[*]};
do
    OutFileName=${WorkDir}/${bam:2:-11}
    echo "bamCoverage -b ${bam} -o ${OutFileName}.bw"
    bamCoverage -b ${bam} -o ${OutFileName}.bw
    
    if [ ! -f ${OutFileName}_count.cntTable ];then
        echo ${OutFileName}
        echo "TEcount --sortByPos \
            --format BAM \
            --mode multi \
            -b ${bam} \
            --project ${OutFileName}_count \
            --GTF ${GENE_GTFFILE} \
            --verbose 0 \
            --TE ${TE_GTFFILE} > ${OutFileName}_count.log"
            
        TEcount --sortByPos \
            --format BAM \
            --mode multi \
            -b ${bam} \
            --project ${OutFileName}_count \
            --GTF ${GENE_GTFFILE} \
            --verbose 0 \
            --TE ${TE_GTFFILE} > ${OutFileName}_count.log  &
        i=$(expr ${i} + 1)
    fi
#    if ( expr ${i} % 10  == 9 )
#    then
        #break
#        sleep 100s
#    fi
done

end_date=`date`
echo "TEcount Start at ${start_date}" #| mail -s "Project Finished at ${end_date}" lux@gwu.edu

