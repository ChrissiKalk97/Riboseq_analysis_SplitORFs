# Riboseq validation of SplitORFs

The SplitORF unique regions are validated with reads from Riboseq data. The Riboseq reads are therefore aligned against the transcriptome or the genome and then intersected with the genomic or transcriptomic coordiantes
of the unique regions.

---

## Directory Structure

This section describes the main directories in the project and what they contain.

```plaintext
Riboseq
│
├── Input_scripts/
├── Uniqueness_scripts/
├── deduplicate/
├── Input_scripts/
├── genomic_alignment/
│   ├── analyze_data.py
│   ├── plot_results.R
│   └── additional_info.txt
├── old_scripts/
└── play_around_bowtie_option_scripts/
