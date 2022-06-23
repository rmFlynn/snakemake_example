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

cd /home/projects-wrighton-2/projects-flynn/may_24_22_meta_g_pipline/RMNP_pipline

snakemake --profile slurm -j 4 -c 20 --dag | dot -Tpdf > dag.pdf
snakemake --profile slurm -j 4 -c 30 --rerun-incomplete
