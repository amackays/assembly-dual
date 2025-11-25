#!/bin/sh
#SBATCH --output=assembly-test.%j.out
#SBATCH --partition=scavenger
#SBATCH --time=4-00:00:00

#source /hpc/home/apm58/.zshrc
source activate /hpc/group/wraylab/apm58/miniconda3/envs/snakemake

snakemake --workflow-profile slurm_general --rerun-incomplete
