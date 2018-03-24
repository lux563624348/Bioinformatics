Output="/home/lxiang/cloud_research/PengGroup/XLi/Data/Paul/34bc/MM10_pMXs_combo"
MM10_FA="/home/lxiang/cloud_research/PengGroup/XLi/Annotation/MM10/mm10.fa"
pMXs_Klf4_FA="/home/lxiang/Raw_Data/Paul/34bc/Reference_gene_seq/pMXs/pMXs-Klf4-full.fa"
pMXs_Oct4_FA="/home/lxiang/Raw_Data/Paul/34bc/Reference_gene_seq/pMXs/pMXs-Oct4-full.fa"
pMXs_Sox2_FA="/home/lxiang/Raw_Data/Paul/34bc/Reference_gene_seq/pMXs/pMXs-Sox2-full.fa"

echo "MM10_pMXs_Combo Build from Bowtie2"

echo "bowtie2-build -f ${MM10_FA},${pMXs_Klf4_FA},${pMXs_Oct4_FA},${pMXs_Sox2_FA}, mm10_pMXs_combo"
bowtie2-build -f ${MM10_FA},${pMXs_Klf4_FA},${pMXs_Oct4_FA},${pMXs_Sox2_FA}, mm10_pMXs_combo
