#!/bin/bash

EXEDIR=/opt/tools

echo "java -Xmx16g -jar $EXEDIR/juicer_tools_1.11.09_jcuda.0.8.jar pre -d Tcf1_KO_na_CD8_Juicebox_input.txt Tcf1_KO_na_CD8_Juicebox.hic mm9"
java -Xmx16g -jar $EXEDIR/juicer_tools_1.11.09_jcuda.0.8.jar pre -d Tcf1_KO_na_CD8_Juicebox_input.txt Tcf1_KO_na_CD8_Juicebox.hic mm9

echo "done"
