#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load fastp
OUTDIR="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples"
heart_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/heart_data/fastq"
leukemia_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq"
endothel_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/bam"

samples=($heart_dir/"ERR3367798.fastq" $heart_dir/"ERR3367797.fastq")

for i in {1..6}; do
  samples+=($endothel_dir"/OHMX20220060_00$i.fastq")
done 

for i in {8..9}; do
  samples+=($leukemia_dir"/SRR1129460$i.fastq")
done

for i in {10..11}; do
  samples+=($leukemia_dir"/SRR112946$i.fastq")
done


for i in "${samples[@]}"
do
SAMPLE=$(basename "$i" .fastq)
FQ=$i
fastp \
    -i ${FQ} \
    -o ${OUTDIR}/${SAMPLE}_fastp.fastq \
    --json ${OUTDIR}/${SAMPLE}.fastp.json \
    --thread 16 \
    --length_required 20\
    --length_limit 45
done
