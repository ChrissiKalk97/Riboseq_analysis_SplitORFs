#!/bin/bash
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate Riboseq

reference_index="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Output/Bowtie_contaminant_index"
output_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Output/Bowtie_contaminant_mapping"
heart_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/ERR3367797_8_fastp.fastq"
leukemia_control_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/leukemia_control_10_11_fastp.fastq"
leukemia_KD_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/leukemia_KD_08_09_fastp.fastq"
threads=32
reference_fasta="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/contaminant_r_t_sn_RNA.fa"

if [ ! -d $reference_index ];then
	mkdir $reference_index
fi

if [ ! -d $output_dir ];then
        mkdir $output_dir
fi

bowtie2-build --threads $threads $reference_fasta $reference_index

samples=("$heart_ribo" "$leukemia_control_ribo" "$leukemia_KD_ribo")


for i in "${samples[@]}"
do
bamfile=${output_dir}/$( basename "$i" .fastq)_contaminants.bam
bowtie2 --threads $threads -x $reference_index -U $i | samtools view -@ $threads -bS > $bamfile
done

