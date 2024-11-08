#!/bin/bash
export PATH=$PATH:/home/fuchs/agschulz/kalk/scripts/nanopore_ana/short_reads/sratoolkit.3.0.7-ubuntu64/bin
#fasterq-dump SRR11294611 --outdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq
#fasterq-dump SRR11294610 --outdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq
fasterq-dump ERR3367798 --outdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/heart_data/fastq
fasterq-dump ERR3367797 --outdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/heart_data/fastq
