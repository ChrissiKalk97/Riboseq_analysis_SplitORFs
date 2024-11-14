#!/bin/bash
#SBATCH --job-name=canriboseq
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=23:40:00
#SBATCH --mail-type=FAIL
#SBATCH --account=agschulz
#SBATCH --mem=120Gb
srun transcriptomic_align_3_primes_canonical.sh
