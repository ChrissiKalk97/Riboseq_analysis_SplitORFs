#!/bin/bash
#SBATCH --job-name=3primes
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=08:40:00
#SBATCH --mail-type=FAIL
#SBATCH --account=agschulz
#SBATCH --mem=120Gb
srun get_3prime_coords_filter_prot_cod_and_tsl.sh

