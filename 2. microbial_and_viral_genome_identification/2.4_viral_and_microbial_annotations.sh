#!/bin/bash
#SBATCH --nodes=1                  # Number of requested nodes
#SBATCH --ntasks=10               # Cores per node requested
#SBATCH --time=7-50:00:00               # Timelimit
#SBATCH --mem=170gb                # Job memory (1024gb=1tb)
#SBATCH --mail-type=BEGIN,END,FAIL         
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

#Viral genome annotations
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.0

#annotate
DRAM-v.py annotate -i 1230_vMAG_scaffolds.fasta -v virsorter2_viral-affi-contigs-for-dramv.tab \
        -o ./output_uniref_1230 --min_contig_size 0 --threads 20  --use_uniref &> log_erpe_dramv_1230_viruses_uniref.txt

#Microbial genome annotations
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.0

#annotate
DRAM.py annotate -i ./MAG_directory/"*.fa" -o output_DRAM --min_contig_size 2500 --threads 25 --use_uniref --use_camper --use_fegenie --use_sulphur --use_vogdb &> log_erpe_dram_1.4.0.txt

DRAM.py distill -i output_DRAM/annotations.tsv -o output_DRAM/distilled