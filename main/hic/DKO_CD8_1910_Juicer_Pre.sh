#!/bin/bash

EXEDIR=/opt/tools

echo "java -Xmx16g -jar $EXEDIR/juicer_tools_1.11.09_jcuda.0.8.jar pre -d DKO_CD8_1910_Juicebox_input.txt.gz DKO_CD8_1910_Juicebox.hic mm9"
java -Xmx16g -jar $EXEDIR/juicer_tools_1.11.09_jcuda.0.8.jar pre -d DKO_CD8_1910_Juicebox_input.txt.gz DKO_CD8_1910_Juicebox.hic mm9

echo "done"
