#!/bin/bash
#SBATCH --nodes=1                  # Number of requested nodes
#SBATCH --ntasks=25               # Cores per node requested
#SBATCH --time=7-50:00:00        	# Timelimit
#SBATCH --mem=140gb                # Job memory (1024gb=1tb)
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

#Run phist
phist.py -k 25 -t 25 ./virus ./host common_kmers.csv predictions.csv

#Run VirHostMatcher
python /opt/VirHostMatcher/vhm.py -v virus -b host -o output

#Run CRASS
#Instructions: Simply run using a bash call in the directory that CRASS output can be found. 

ls *spacers.gv | awk -F "_" '{ print $2}' > splitfiles.txt #Make gv files into a line separated list.

mkdir crass_pipeliner
cd crass_pipeliner

split -l 20 ../splitfiles.txt split_ #Split comma separated list into individual files of 20 spacers each.

ls split_* > split_file_paths.txt #Make the directory for final sed commands.

while read x; do
	sed '$!s/$/,/' "$x" > "$x"_comma
	sed '{:q;N;s/\n//g;t q}' "$x"_comma > "$x"_comma_list
done <split_file_paths.txt

cat *_comma_list > comma_list_concatenated.txt #make a concatenated file of each of the split files where each file content is a single line. Then we feed each of these lines as the argument for the crisprtools stat command.

#run the crisprtool stats
while read j; do
		crisprtools stat -g "$j" -apH ../crass.crispr > crisprtools_"$j"_stats_sampleID.out
done <comma_list_concatenated.txt

#now run the crisprtools extract
while read j; do
		crisprtools extract -o "crisprtools_extract" -g "$j" -s -xC -d -f ../crass.crispr
		cat crisprtools_extract/*_direct_repeats.fa > All_direct_repeats.fa
		cat crisprtools_extract/*_spacers.fa > All_spacers.fa
done <comma_list_concatenated.txt

#Now that we have our crispr direct repeats parsed, we need to blast them to our MAGs. I have the MAG database of 125 MAGs with unique header ID's here: /home/projects-wrighton-2/GROWdb/WHONDRS_2018/metaG/2018_rivers/erpe/erpe_JRR/surface_and_pore_viral_db_analyses/crass/125_erpe_MAGs_no2019_ge5kb.fa. We want to do this only on scaffolds that are greater than 5kb. We went from 35443 to 13709 by using the 5kb cutoff.