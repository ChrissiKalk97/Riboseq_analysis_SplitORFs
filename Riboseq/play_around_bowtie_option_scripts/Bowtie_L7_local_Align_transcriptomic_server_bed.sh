#!/bin/bash

#Help message:
usage="
Usage: ./Bowtie_Align_transcriptomic.sh [-options] numberOfThreads BowtieBaseName Reads.fastq out unique_regions.bed exonAnnotation genomicAnnotation

numberOfThreads 		Int setting the number of threads to use with BOWTIE2
BowtieBaseName 			BaseName to use for creating BOWTIE2 Index files, or if allready created, BaseName of the existing index files
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz)
out					Base name of the output files
unique_regions.bed		BED-file containing the unique_DNA_regions up for validation
transcripts.fa	fasta file with reference transcripts (NMD, RI etc.)

where:
-h			show this help
-i transcripts.fa	create new index files for the provided transcripts.fa"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':hi' option; do
  case "$option" in
    h) echo "$usage"
       exit 1
       ;;
	i) echo "creating new index files"
	   if [[ $# -ne 8 ]]; then #check for right number of arguments
		echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
		echo "$usage" >&2
		exit 1
	   fi
       bowtie2-build --threads $3 $2 $4
	   ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
  shift 2
done

if [[ $# -ne 6 ]]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
fi

numberOfThreads=$1
bowtieBaseName=$2
Riboreads=$3
bamfile=$4.bam
unique_regions=$5
coordinates_3_prime=$6

out_path=$(dirname $bamfile)
#align Riboreads against the transcriptome
bowtie2 --threads $numberOfThreads -L 7 --local -x $bowtieBaseName -U $Riboreads | samtools view -@ $numberOfThreads -bS > $bamfile
bedfile=$out_path/$(basename $bamfile .bam).bed

echo "converting bam to bed"
bedtools bamtobed -i $bamfile > $bedfile
sortedBedfile=$out_path/$(basename $bedfile .bed)_intersect_counts_sorted.bed

echo "intersecting with unique regions"
bedtools intersect -c -F 0.33 -a $unique_regions -b $bedfile | sort -nr -k 4 > $sortedBedfile
intersectbedfilerelativesorted=$out_path/$(basename $bedfile .bed)_intersect_counts_relative_sorted.bed
cat $sortedBedfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$5, $5/($3-$2)}' | sort -nr -k 6 > $intersectbedfilerelativesorted 

echo "sorting bamfile and indexing"
sortedbamfile=$(echo $bamfile | rev | cut -f 2- -d '.' | rev)_sorted.bam
samtools sort -@ 10 $bamfile > $sortedbamfile
samtools index -@ 10 $sortedbamfile


echo "Calculating random regions from 3 prime UTRs"
randomfile=$out_path/$(basename $bamfile .bam)_random_background_regions.bed
python ./Uniqueness_scripts/BackgroundRegions_bed.py $unique_regions $coordinates_3_prime $randomfile
randomintersectfile=$out_path/$(basename $bamfile .bam)_random_intersect_counts.bed
bedtools intersect -c -F 0.33 -a $randomfile -b $bedfile > $randomintersectfile
randomintersectfilesorted=$out_path/$(basename $bamfile .bam)_random_intersect_counts_relative_sorted.bed
cat $randomintersectfile | awk -v OFS='\t' '{print $1,$2,$3,$4, $4/($3-$2)}' | sort -n -r -k 5 > $randomintersectfilesorted
