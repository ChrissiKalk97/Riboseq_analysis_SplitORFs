# Riboseq validation of SplitORFs

The SplitORF unique regions are validated with reads from Riboseq data. The Riboseq reads are therefore aligned against the transcriptome or the genome and then intersected with the genomic or transcriptomic coordiantes
of the unique regions.

---

## Directory Structure

The main directory contains subdectories with collections of scripts as well as the main scripts (latest) used for the transcriptomic alignment and evaluation of the transcriptomic validation method.<br>
The sbatch_transcriptomic_align_08_10_24_new_3_primes_new_unique_regions.sh needs to be called in order to run the transcriptomic analysis.<br>
The Input_scripts folder contains scripts from prefetching the Input_data until their pre-processing with fastp.<br>
The 

```plaintext
Riboseq
│
├── Input_scripts/
├── Random_regions/
├── deduplicate/
├── Input_scripts/
├── genomic_alignment/
│   ├── analyze_data.py
│   ├── plot_results.R
│   └── additional_info.txt
├── old_scripts/
├── play_around_bowtie_option_scripts/
├── Bowtie_Align_transcriptomic_server_bed.sh
├── RiboSeqReportTranscriptomic_empirical_dist.Rmd
├── helper_functions_Riboseq_report_new.R
├── sbatch_transcriptomic_align_08_10_24_new_3_primes_new_unique_regions.sh
└── transcriptomic_align_08_10_24_new_3_primes_new_unique_regions.sh

