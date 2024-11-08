#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load fastqc
if [ ! -d /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples/fastqc ]; then
	mkdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples/fastqc
fi
fastqc -o /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples/fastqc -t 20 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples/*.fastq

