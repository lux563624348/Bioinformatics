#\!/bin/bash
# gzip .sam, .bed, .wig files
find $1 \( -name '*.sam' \) -exec gzip --verbose {} \;
find $1 \( -name '*.bed' \) -exec gzip --verbose {} \;
find $1 \( -name '*.wig' \) -exec gzip --verbose {} \;
