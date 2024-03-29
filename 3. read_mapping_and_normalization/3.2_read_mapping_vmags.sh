#!/bin/bash

# pre-processing:
# drep vMAGs at 95% ANI across 85% of the shortest contig.
# concatenate vMAGs into single file for mapping reference.
# list.txt is a list of trimmed read names.
# flags for CoverM and %ID based off of Roux et al., 2017 https://peerj.com/articles/3817/ and Gregory et al., 2019 https://www.sciencedirect.com/science/article/pii/S0092867419303411#sec3
# bowtie2 snippets and reformat.sh gotten from Bridget McGivern after she spent a lot of time figuring out best mapping protocols.

cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2

bowtie2-build /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/1250_vMAGs/all_viral_genomes_erpe_1250.fna vmag_DB --large-index --threads 10 #create bowtie2 index with a concatenated viral genomes fasta.

for element in $(<$1)
do
echo "begin bowtie2 '$element' mapping"
mkdir "$element"_mapping
cd "$element"_mapping
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 15 -x /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/vmag_DB -S "$element".sam -1 /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/trimmed_read_copies_for_mapping/"$element"_R1.fastq -2 /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/trimmed_read_copies_for_mapping/"$element"_R2.fastq #run bowtie2 with indexed genomes. This is basically preset of --fast option specified in bowtie2 manual. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

echo "begin samtools"
samtools view -@ 15 -bS "$element".sam > "$element".bam #sam to bam
reformat.sh -Xmx100g idfilter=0.95 pairedonly=t primaryonly=t in="$element".bam out="$element"_mapped_95ID.FILTERED.bam #this will then give you the actual 95% ID mapping as well as remove the redundancies of bowtie2 mapping.
samtools sort -T "$element".95ID.sorted -o "$element".95ID.sorted.bam "$element"_mapped_95ID.FILTERED.bam -@ 15 #sort the bam file
mv "$element".95ID.sorted.bam /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1230_final_set_nov1_22_bowtie2/sorted_BAM_95ID #now we move out the sorted BAM into its own directory where all will be for coverM. 
echo "Finished mapping and samtools for '$element' in list" 
cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2 #move back up to root and restart the process.
done

#now that we've moved over all the sorted SAMs into their own directory, we run coverM.

cd /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/sorted_BAM_95ID

echo "begin coverM"
coverm genome -m mean --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/1250_vMAGs --bam-files *.sorted.bam --min-covered-fraction 75 -t 15 --output-format dense &> coverM_min-covered75.txt #get 75% minimum covered fraction
coverm genome -m reads_per_base --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/1250_vMAGs --bam-files *.sorted.bam -t 15 --min-covered-fraction 0 --output-format dense &> coverM_reads.per.base.txt #to get coverage / depth (need to multiply resulting value by 151bp)
coverm genome -m count --genome-fasta-extension fasta --genome-fasta-directory /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/read_mapping/vMAGs/1250_final_set_010422_bowtie2/1250_vMAGs --bam-files *.sorted.bam --min-covered-fraction 0 -t 15 --output-format dense &> coverM_counts.txt #Get counts for geTMM
echo "coverM coverage, depth, and counts finished here"
echo "finished mapping scripts"
