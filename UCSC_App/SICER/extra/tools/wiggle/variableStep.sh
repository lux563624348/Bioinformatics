#!/bin/sh
##RNAtail seq. peaks are sharp. thus span region is only 10 bps. 
SPAN=10
S=$1
R=$2
rm -f ${R}
head -45 ${S} | egrep "^browser|^track" > ${R}
echo "track type=wiggle_0 name=$3" >> ${R}
grep "^chr" ${S} | cut -f1 | sort -u > chr.list
cat chr.list | while read C
do
    echo "variableStep chrom=${C} span=${SPAN}" >> ${R}
    awk '{if (match($1,"^'"${C}"'$")) { print } }' ${S} | sort -k2n | awk '
{
    printf "%d\t%g\n", $2+1, $4
}
' >> ${R}
done
