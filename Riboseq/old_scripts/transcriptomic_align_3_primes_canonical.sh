#!/bin/bash -l
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
nmd="./Output/Bowtie_transcriptome_heart_leuk_Siragusa_canonical/NMD_transcriptome"
unique_region_dir="/scratch/fuchs/agschulz/kalk/k4neo/unique_regions"
ri="./Output/Bowtie_transcriptome_heart_leuk_Siragusa_canonical/RI_transcriptome"

outputBowtie="./Output/Bowtie_transcriptome_heart_leuk_Siragusa_canonical"
#The following are the paths to the riboseq reads
#the two samples are merged
heart_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/heart_data/fastq/ERR3367797_cropped_38.fastq"
leukemia_control_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq/leukemia_control_10_11.fastq"
leukemia_KD_ribo="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/leukemia_data/fastq/leukemia_KD_08_09.fastq"
#controls
control="./bam/control.fastq"
#treatments
treatment="./bam/treatment_wo_05.fastq"

#file with the 3'UTR background regiosn
three_primes="./Input/3_primes_Ens_110_canonical_tsl.bed"

#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputBowtie/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against transcripts"
source ./Bowtie_Align_transcriptomic_server_bed.sh -i ./Input/cDNA_fasta_Ensembl_transcripts_110.fa 10 $outputBowtie/Transcriptomic_Bowtie_index $heart_ribo $nmd/heart_NMD\
 $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $heart_ribo\
 $ri/heart_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes
echo "===================	heart data finished"


source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_control_ribo\
 $nmd/leukemia_control_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_control_ribo\
 $ri/leukemia_control_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes

echo "===================       leukemia control data finished"

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_KD_ribo\
 $nmd/leukemia_KD_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $leukemia_KD_ribo\
 $ri/leukemia_KD_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes

echo "===================       leukemia KD data finished"

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
 $nmd/endothel_control_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
 $ri/endothel_control_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes


source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $treatment\
 $nmd/endothel_treatment_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

source ./Bowtie_Align_transcriptomic_server_bed.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $treatment\
 $ri/endothel_treatment_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes
echo "===================       Endothel data inished"

