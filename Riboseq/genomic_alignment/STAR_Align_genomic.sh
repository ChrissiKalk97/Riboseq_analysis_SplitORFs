#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load star@2.7.10b

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
StarIndex=$2
Riboreads=$3
bamfile=${4}.bam
unique_regions=$5
coordinates_3_prime=$6

out_path=$(dirname $bamfile)
#align Riboreads against the genome
STAR\
 --runThreadN $numberOfThreads\
 --alignEndsType EndToEnd\
 --outSAMstrandField intronMotif\
 --alignIntronMin 20\
 --alignIntronMax 1000000\
 --genomeDir $StarIndex\
 --readFilesIn $Riboreads\
 --twopassMode Basic\
 --outFileNamePrefix $out_path/$(basename $bamfile .bam)
#/scratch/fuchs/agschulz/kalk/star/reference_110_ribo
samtools view -@ $numberOfThreads -bo $bamfile $out_path/$(basename $bamfile .bam)Aligned.out.sam
sortedBamFile=$out_path/$(basename $bamfile .bam)_sorted.bam
samtools sort -o $sortedBamFile $bamfile
samtools index -@ 10 $sortedBamFile
bedfile=$out_path/$(basename $bamfile .bam).bed

echo "converting bam to bed"
bedtools bamtobed -bed12 -i $sortedBamFile > $bedfile


echo "intersecting with unique regions"
sortedBedfile=$out_path/$(basename $bedfile .bed)_intersect_counts_sorted.bed
sorted_unique_regions=$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed

# sort unique regions for intersect
sort -k1,1 -k2,2n $unique_regions > $sorted_unique_regions

# check if genome file for chromosome order for sorted intersect exists
# if not create
if [ ! -e  $out_path/genome_file.txt ]; then
    cut -f1,2 /scratch/fuchs/agschulz/kalk/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa.fai >\
    $outpath/genome_chrom_ordering.txt
fi

# intersect such that both entries are reported
# entry A is always reported with a null B feature, if no intersect
# file A the unique regions is a 6-bed, the Riboreads are a 12-bed
# the last col 19 gives the number of basepairs of overlap

bedtools intersect\
  -s\
  -wao\
  -F 0.33\
  -a $sorted_unique_regions\
  -b $bedfile\
  -sorted\
  -g $outpath/genome_chrom_ordering.txt\
  > $sortedBedfile\
  | sort -nr -k19


# intersectbedfilerelativesorted=$out_path/$(basename $bedfile .bed)_intersect_counts_relative_sorted.bed
# cat $sortedBedfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$5, $5/($3-$2)}' | sort -nr -k 6 > $intersectbedfilerelativesorted 


echo "Calculating random regions from 3 prime UTRs"
randomfile=$out_path/$(basename $bamfile .bam)_random_background_regions.bed
python /home/fuchs/agschulz/kalk/scripts/SplitORFs/Riboseq/Uniqueness_scripts/BackgroundRegions_bed_genomic.py\
 $unique_regions\
 $coordinates_3_prime\
 $randomfile

randomintersectfile=$out_path/$(basename $bamfile .bam)_random_intersect_counts.bed
bedtools intersect -s -wao -F 0.33 -a $randomfile -b $bedfile > $randomintersectfile
# randomintersectfilesorted=$out_path/$(basename $bamfile .bam)_random_intersect_counts_relative_sorted.bed
# cat $randomintersectfile | awk -v OFS='\t' '{print $1,$2,$3,$4, $4/($3-$2)}' | sort -n -r -k 5 > $randomintersectfilesorted
