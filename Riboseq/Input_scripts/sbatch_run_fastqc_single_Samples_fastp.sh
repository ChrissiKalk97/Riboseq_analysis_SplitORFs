#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=23:40:00
#SBATCH --mail-type=FAIL
#SBATCH --account=agschulz
#SBATCH --mem=120gb
srun run_fastqc_single_Samples_fastp.sh
