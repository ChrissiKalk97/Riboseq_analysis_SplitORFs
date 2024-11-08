#!/bin/bash -l
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
outputBowtie="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Output/Bowtie_transcriptome_08_10_24_new_unique_new_3_prime"
nmd=${outputBowtie}"/NMD_transcriptome"
unique_region_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/Unique_regions/all_trans_refseq_tsl_filtered_ref"
ri=${outputBowtie}"/RI_transcriptome"

if [ ! -d $outputBowtie ];then
	mkdir $outputBowtie
fi

if [ ! -d $nmd ];then
        mkdir $nmd
fi
if [ ! -d $ri ];then
        mkdir $ri
fi


#The following are the paths to the riboseq reads
#the two samples are merged
heart_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/ERR3367797_8_fastp.fastq"
leukemia_control_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/leukemia_control_10_11_fastp.fastq"
leukemia_KD_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/leukemia_KD_08_09_fastp.fastq"
#controls
control="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/control_fastp.fastq"
#treatments
treatment="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/treatment_wo_05_fastp.fastq"

#file with the 3'UTR background regiosn
three_primes="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/three_primes_tsl1_2_refseq_prot_cod.bed"

#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputBowtie/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against transcripts"
#source ./Bowtie_Align_transcriptomic_server_bed.sh -i /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/cDNA_fasta_Ensembl_transcripts_110.fa 10 $outputBowtie/Transcriptomic_Bowtie_index $heart_ribo $nmd/heart_NMD\
# $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $heart_ribo\
# $ri/heart_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes
echo "===================	heart data finished"


#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_control_ribo\
# $nmd/leukemia_control_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_control_ribo\
# $ri/leukemia_control_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes

echo "===================       leukemia control data finished"

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_KD_ribo\
# $nmd/leukemia_KD_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_KD_ribo\
# $ri/leukemia_KD_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes

echo "===================       leukemia KD data finished"

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
# $nmd/endothel_control_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

#source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
# $ri/endothel_control_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes


source ./Bowtie_Align_transcriptomic_server_bed.sh  -i /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/cDNA_fasta_Ensembl_transcripts_110.fa 10 $outputBowtie/Transcriptomic_Bowtie_index $treatment\
 $nmd/endothel_treatment_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $treatment\
 $ri/endothel_treatment_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes
echo "===================       Endothel data inished"

