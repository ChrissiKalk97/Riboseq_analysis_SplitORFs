#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load fastp
OUTDIR=/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp
heart_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/heart_data/fastq"
leukemia_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq"
endothel_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/bam"

samples=($heart_dir/"ERR3367797_8.fastq" $leukemia_dir"/leukemia_control_10_11.fastq" $leukemia_dir"/leukemia_KD_08_09.fastq"\
 $endothel_dir"/control.fastq" $endothel_dir"/treatment_wo_05.fastq")


for i in "${samples[@]}"
do
SAMPLE=$(basename "$i" .fastq)
FQ=$i
fastp \
    -i ${FQ} \
    -o ${OUTDIR}/${SAMPLE}_fastp.fastq \
    --json ${OUTDIR}/${SAMPLE}.fastp.json \
    --thread 32 \
    --length_required 20
done
