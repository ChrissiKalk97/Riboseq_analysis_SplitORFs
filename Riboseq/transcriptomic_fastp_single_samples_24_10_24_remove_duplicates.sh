#!/bin/bash -l
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
outputBowtie="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Output/Bowtie_transcriptome_24_10_24_single_samples_dedup"
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
input_data="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/fastp/fastp_single_samples/*.fastq"

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


# bowtie2-build --threads 16 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/cDNA_fasta_Ensembl_transcripts_110.fa $outputBowtie/Transcriptomic_Bowtie_index

# for i in $input_data; do
# sample_name=$(basename $i _fastp.fastq)
# source ./Bowtie_Align_transcriptomic_remove_duplicates.sh  16 $outputBowtie/Transcriptomic_Bowtie_index $i\
#  $nmd/${sample_name}_NMD $unique_region_dir/Unique_DNA_Regions_for_riboseq_NMD.bed $three_primes

# source ./Bowtie_Align_transcriptomic_remove_duplicates.sh 16 $outputBowtie/Transcriptomic_Bowtie_index $i\
#  $ri/${sample_name}_RI $unique_region_dir/Unique_DNA_Regions_for_riboseq_RI.bed $three_primes

# echo "===================       Sample $sample_name finished"

# done


source activate my_r_env


Rscript -e "if (!requireNamespace('rmarkdown', quietly = TRUE)) install.packages('rmarkdown', repos='http://cran.us.r-project.org')"

R -e "library(rmarkdown); rmarkdown::render(input = 'RiboSeqReportTranscriptomic_empirical_dist_dedup.Rmd', output_file = '$outputBowtie/Riboseq_report.pdf', params=list(args = c('$outputBowtie')))"

