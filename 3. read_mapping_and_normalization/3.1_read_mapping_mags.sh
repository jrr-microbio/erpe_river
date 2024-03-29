#!/bin/bash

# pre-processing:
# concatenate MAGs into single file for mapping reference.
#Select medium-high quality genomes
# list.txt is a list of trimmed read names.

cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs

bowtie2-build /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/erpe_2018_drepped_125_MAGs.fasta mag_DB --large-index --threads 10 #create bowtie2 index with a concatenated genomes fasta.

for element in $(<$1)
do
echo "begin bowtie2 '$element' mapping"
mkdir "$element"_mapping
cd "$element"_mapping
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 15 -x /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/mag_DB -S "$element".sam -1 /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/trimmed_read_copies_for_mapping/"$element"_R1.fastq -2 /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/trimmed_read_copies_for_mapping/"$element"_R2.fastq #run bowtie2 with indexed genomes. This is basically preset of --fast option specified in bowtie2 manual. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

echo "begin samtools"
samtools view -@ 15 -bS "$element".sam > "$element".bam #sam to bam
reformat.sh -Xmx100g idfilter=0.95 pairedonly=t primaryonly=t in="$element".bam out="$element"_mapped_95ID.FILTERED.bam #this will then give you the actual 95% ID mapping. This was wrong before on bbmap as the 95% ID cutoffs were not actually cutoffs. UPDATE: Bowtie2 is "leaky" in the sense that when a read mapping hit is equally good, it will report all hits. This additional flag ensures that you are only using non-ambiguous calls and truly competitive mapping as well as making sure that all your mapping outputs are from actually paired reads that mapped and not single/unpaired sets.
samtools sort -T "$element".95ID.sorted -o "$element".95ID.sorted.bam "$element"_mapped_95ID.FILTERED.bam -@ 15 #sort the bam file
mv "$element".95ID.sorted.bam /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/sorted_BAM_95ID #now we move out the sorted BAM into its own directory where all will be for coverM. 
echo "Finished mapping and samtools for '$element' in list" 
cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs1250_final_set_010422_bowtie2 #move back up to root and restart the process.
done

#now that we've moved over all the sorted SAMs into their own directory, we run coverM.

cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/sorted_BAM_95ID

echo "begin coverM"
coverm genome -m mean --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/dereplicated_genomes --bam-files *.sorted.bam --min-covered-fraction 75 -t 15 --output-format dense &> coverM_min-covered75.txt #get 75% minimum covered fraction
coverm genome -m reads_per_base --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/dereplicated_genomes --bam-files *.sorted.bam -t 15 --min-covered-fraction 0 --output-format dense &> coverM_reads.per.base.txt #to get coverage / depth (need to multiply resulting value by 151bp)
coverm genome -m count --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/MAGs/dereplicated_genomes --bam-files *.sorted.bam --min-covered-fraction 0 -t 15 --output-format dense &> coverM_counts.txt #Get counts for geTMM
echo "coverM coverage, depth, and counts finished here"
echo "finished mapping scripts"

###################################New addition
#Bowtie2 mapping is competitive by default - but is also leaky. When it gets 2 identical hits, it reports all hits and so you get hits that are redundant. In order to remove this, we apply the additional reformat.sh command to the 95ID filtered and SORTED.bam files for each set of reads:
reformat.sh -Xmx100g pairedonly=t primaryonly=t in="$element".95ID.sorted.bam out="$element"_mapped_95ID.FILTERED_sorted_comp-map.bam
#Rerun coverM on the output file.
