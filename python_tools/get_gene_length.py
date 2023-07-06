#!/usr/bin/env python3
import re
import sys

wantedGenes = []
with open(sys.argv[1], mode='r') as f:
    wantedGenes = f.readlines()

wantedGenes = [ x.rstrip() for x in wantedGenes]

print('# Gene\tstart\tend\tlength')
with open(sys.argv[2], mode='r') as f:
    for line in f:
        if not line.startswith('#'): 
            fields = line.split('\t')
            annotations = fields[8].split(' ')
            geneName = re.sub('[";]', '', annotations[1])
            if geneName in wantedGenes and fields[2] == 'gene':
                print('%s\t%s\t%s\t%d' %
                      (geneName, fields[3],fields[4],
                       int( int(fields[4]) - int(fields[3]) )
                      ))
