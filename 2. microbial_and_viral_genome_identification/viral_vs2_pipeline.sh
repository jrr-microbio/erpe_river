#!/bin/bash

#$1 is a list of directories that contains assembly names

###########VirSorter2, checkv and virsorter2 2nd run

source /opt/Miniconda2/miniconda2/bin/activate virsorter2
cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/pore_water/

for element in $(<$1)
do
mkdir "$element".virsorter2out
cd "$element".virsorter2out
virsorter run -w vs2_"$element"_pass1 --keep-original-seq -i ../"$element"/"$element".contigs.fa -w vs2-pass1 --include-groups dsDNAphage,ssDNA --min-length 10000 --min-score 0.5 -j 45 all
checkv contamination /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/pore_water/"$element".virsorter2out/vs2-pass1/final-viral-combined.fa checkv_out
checkv completeness /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/pore_water/"$element".virsorter2out/vs2-pass1/final-viral-combined.fa checkv_out
checkv repeats /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/pore_water/"$element".virsorter2out/vs2-pass1/final-viral-combined.fa checkv_out
checkv quality_summary /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/pore_water/"$element".virsorter2out/vs2-pass1/final-viral-combined.fa checkv_out
virsorter run --seqname-suffix-off --viral-gene-enrich-off --prep-for-dramv -i checkv_out/cleaned_contigs.fna -w vs2-pass2 --include-groups dsDNAphage,ssDNA --min-length 10000 --min-score 0.5 -j 45 all
cd ../
done
echo "VirSorter2 and checkV are complete"

######Once VirSorter2 is complete on all assemblies, concatenate the genomes from all and cluster as below.

Cluster_genomes_5.1.pl -f /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/all_erpe_vMAGs_final-viral-combined.fa -c 85 -i 95 -t 6
