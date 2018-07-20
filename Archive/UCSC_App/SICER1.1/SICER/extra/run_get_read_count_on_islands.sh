#!/bin/bash
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng 
#

SICER=/home/data/SICER1.1/SICER

echo "sh $SICER/extra/get_read_count_on_islands.sh hg18 150 72h_EM_EZh2-W200-G600-E500-islandfiltered.bed 72h_EM_EZh2-W200-G600-E500.scoreisland 72h_EM_EZh2-W200-G600-E500.islands"

sh $SICER/extra/get_read_count_on_islands.sh hg18 150 72h_EM_EZh2-W200-G600-E500-islandfiltered.bed 72h_EM_EZh2-W200-G600-E500.scoreisland 72h_EM_EZh2-W200-G600-E500.islands . .