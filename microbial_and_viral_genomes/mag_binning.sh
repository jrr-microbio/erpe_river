#!/bin/bash

#$1 list of names of the assemblies.

#######mapping reads to scaffolds greater than 2500bp

for element in $(<$1)
do
cd "$element"
bbmap.sh -Xmx48G threads=15 minid=85 overwrite=t ref="$element"_assembly_2500.fa in1=/home/projects/WHONDRS/metaG/2018_rivers/erpe/erpe_JRR/surface_water/trimmed_"$element"_R1.fastq in2=/home/projects/WHONDRS/metaG/2018_rivers/erpe/erpe_JRR/surface_water/trimmed_"$element"_R2.fastq out="$element"_mapped.sam
cd ..
done
echo "Mapping finished"

##########samtools for sorting mapping output file and then binning

for element in $(<$1)
do
cd "$element"
samtools view -@ 15 -bS "$element"_mapped.sam > "$element"_mapped.bam
samtools sort -T "$element"_t_2500.sorted -o "$element"_2500.sorted.bam "$element"_mapped.bam -@ 15
runMetaBat.sh "$element"_assembly_2500.fa "$element"_2500.sorted.bam
cd ..
done
echo "Sorting and Mapping Completed"
