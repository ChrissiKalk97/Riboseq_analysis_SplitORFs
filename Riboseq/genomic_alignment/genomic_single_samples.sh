#!/bin/bash -l
. ~/spack/share/spack/setup-env.sh
spack load star@2.7.10b

PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
outputBowtie="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Output/Bowtie_genomic_16_12_24"
nmd=${outputBowtie}"/NMD_genome"
unique_region_dir="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/Unique_regions/genomic"
ri=${outputBowtie}"/RI_genome"

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
three_primes="/scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/genomic/3_primes_genomic.bed"

#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputBowtie/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against genome"


# STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /scratch/fuchs/agschulz/kalk/star/reference_110_ribo --genomeFastaFiles /scratch/fuchs/agschulz/kalk/star/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa\
#  --sjdbGTFfile /scratch/fuchs/agschulz/kalk/star/Homo_sapiens.GRCh38.110.chr.gtf --sjdbOverhang 28




# for i in $input_data; do
# sample_name=$(basename $i _fastp.fastq)
# source ./STAR_Align_genomic.sh  16 /scratch/fuchs/agschulz/kalk/star/reference_110_ribo $i\
#  $nmd/${sample_name}_NMD $unique_region_dir/Unique_DNA_regions_genomic_NMD_16_12_24.bed $three_primes

# source ./STAR_Align_genomic.sh 16 /scratch/fuchs/agschulz/kalk/star/reference_110_ribo $i\
#  $ri/${sample_name}_RI $unique_region_dir/Unique_DNA_regions_genomic_RI_16_12_24.bed $three_primes

# echo "===================       Sample $sample_name finished"

# done


source activate my_r_env


Rscript -e "if (!requireNamespace('rmarkdown', quietly = TRUE)) install.packages('rmarkdown', repos='http://cran.us.r-project.org')"

R -e "library(rmarkdown); rmarkdown::render(input = 'RiboSeqReportGenomic_empirical_dist.Rmd', output_file = '$outputBowtie/Riboseq_report.pdf', params=list(args = c('$outputBowtie')))"

