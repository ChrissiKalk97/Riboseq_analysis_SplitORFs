#!/bin/bash -l
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate myenvname
python get_3prime_genomic_coords.py\
 /scratch/fuchs/agschulz/kalk/star/filtered_Ensembl_reference/Ensembl_equality_and_TSL_filtered.gtf\
 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/genomic/3_prime_UTR_coords_genomic_Ensembl_110.txt\
 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/genomic/3_primes_genomic.bed

