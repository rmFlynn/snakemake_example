#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=24:00:00
#SBATCH --job-name=snakemake_main
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<your_username>@colostate.edu
#SBATCH --partition=debug
#SBATCH --output=snakemake_main_%j.out


### Release version - Assembly test

# mv ./results/all_counts.tsv ./results/all_counts_first_30.tsv

source /opt/Miniconda2/miniconda2/bin/activate scripts

# Normal run
#> snakemake --profile slurm -j 25 -c 1  --keep-incomplete --notemp
# Debug dry run
snakemake --profile slurm -j 3 -c 30 --keep-incomplete --notemp --dry-run
# Debug run
#> snakemake --profile slurm -j 25 -c 1  --keep-incomplete --notemp
# Make a visual of what you will run
#>snakemake --profile slurm -j 6 -c 20 --rerun-incomplete --keep-incomplete --notemp --dag | dot -Tpdf > dag.pdf
# 
