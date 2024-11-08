#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load fastqc
if [ ! -d /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastqc ]; then
	mkdir /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastqc
fi
fastqc -o /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastqc -t 20 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/*.fastq

