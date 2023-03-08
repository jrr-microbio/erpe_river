#!/bin/bash
#SBATCH --nodes=1                  # Number of requested nodes
#SBATCH --ntasks=15                 # Cores per node requested
#SBATCH --time=7-50:00:00        	# Timelimit
#SBATCH --mem=120gb                # Job memory (1024gb=1tb)
#SBATCH --mail-type=BEGIN,END,FAIL         
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

#Microbial taxonomy using GTDBtk
gtdbtk classify_wf --extension fa --genome_dir ./all_microbial_genomes_all_samples --out_dir ./gtdbtk_all_nondrepped

#Viral taxnomy using vContact2
#NOTE: erpe_river_full_vContact2_db.faa has been deposited in Zenodo due to github size limitations.

source /opt/Miniconda2/miniconda2/bin/activate vContact2

vcontact2_gene2genome -p erpe_river_full_vContact2_db.faa -o ./viral_genomes_g2g.csv -s 'Prodigal-FAA'

vcontact2 -r erpe_river_full_vContact2_db.faa --rel-mode Diamond -p ./viral_genomes_g2g.csv --db ProkaryoticViralRefSeq211-Merged --pcs-mode MCL --vcs-mode ClusterONE --pc-evalue 0.0001 --reported-alignments 25 --max-overlap 0.8 --penalty 2.0 --haircut 0.1 --pc-inflation 2.0 --vc-inflation 2.0 --min-density 0.3 --min-size 2 --vc-overlap 0.9 --vc-penalty 2 --vc-haircut 0.55 --merge-method single --similarity match --seed-method nodes --sig 1.0 --max-sig 300 --mod-inflation 5 --mod-sig 1.0 --mod-shared-min 3 --link-sig 1.0 --link-prop 0.5 --verbose -vv --c1-bin /opt/Miniconda2/miniconda2/envs/vContact2/bin/cluster_one-1.0.jar -o ./output_vContact2 -t 15 

source /opt/Miniconda2/miniconda2/bin/deactivate