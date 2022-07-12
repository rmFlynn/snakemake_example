#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=24:00:00
#SBATCH --job-name=snakemake_main
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mborton@colostate.edu
#SBATCH --partition=debug
#SBATCH --output=snakemake_main_%j.out

source /opt/Miniconda2/miniconda2/bin/activate scripts
snakemake -j 7 -c 1 --rerun-incomplete --notemp --dry-run

# cd /home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/coverm_MetaG_merger_test/mini_pipline

# snakemake --profile slurm -j 4 -c 20 --dag | dot -Tpdf > dag.pdf
snakemake --profile slurm -c 30 --rerun-incomplete --notemp 
snakemake -j 7 -c 1 --rerun-incomplete --notemp


