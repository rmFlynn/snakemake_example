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

cd /home/projects-wrighton-2/projects-flynn/rmnp_pipeline/test_runs_of_pipline/jun_23_metat_check/RMNP_pipline

rm 
snakemake --profile slurm -j 5 -c 20 --rerun-incomplete --keep-incomplete --notemp 
