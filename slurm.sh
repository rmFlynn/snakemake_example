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



# mv ./results/all_counts.tsv ./results/all_counts_first_30.tsv

source /opt/Miniconda2/miniconda2/bin/activate scripts

snakemake --profile slurm -j 20 -c 3  --keep-incomplete --notemp --dry-run
snakemake --profile slurm -j 6 -c 20 --rerun-incomplete --keep-incomplete --notemp --dag | dot -Tpdf > dag.pdf
realpath dag.pdf
snakemake --profile slurm -j 25 -c 1  --keep-incomplete --notemp
# 
