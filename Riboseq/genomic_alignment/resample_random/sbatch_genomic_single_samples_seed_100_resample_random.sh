#!/bin/bash
#SBATCH --job-name=11_20_RiGenResample
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-12:40:00
#SBATCH --mail-type=FAIL
#SBATCH --account=agschulz
#SBATCH --mem=120Gb
srun genomic_single_samples_seed_100_resample_random.sh
