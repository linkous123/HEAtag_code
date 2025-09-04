# HEATag: Code & Pipelines for Cost‑Efficient Sample Labeling in Multiplexed scRNA‑seq

This repository hosts code used in the study  
*“A universal cost‑efficient sample labeling approach for multiplexed single‑cell RNA‑seq based on recombinant HUH‑endonuclease–agglutinin tagging.”*

---

## Table of Contents
- [Overview](#overview)
- [Quick Start](#quick-start)
- [1) HEATag Barcode Counting](#1-heatag-barcode-counting)
- [2) Genome Functional Annotation](#2-genome-functional-annotation)
- [3) Single‑Cell RNA‑seq Processing](#3-single-cell-rna-seq-processing)
- [4) SNP‑based Demultiplexing](#4-snp-based-demultiplexing)
- [5) Visualization](#Visualization)
- [Cite](#cite)
- [License](#license)
- [Contact](#contact)
---

## Overview

This repo provides:
- **HEATag barcode counting** from paired‑end FASTQ files.
- **Functional annotation** via local UniProt/NR BLAST and reciprocal best hits (RBH).
- **Downstream scRNA‑seq analysis** (QC, filtering, demultiplexing, DE & enrichment).
- **SNP‑based demultiplexing** using scSplit/Freemuxlet.
- - **Visualization** visualization code

---

## Quick Start

1. Prepare paired‑end HEATag FASTQs and `barcodes.tsv` (cell barcodes from your count matrix).
2. Edit the **config file** to specify cell/sample‑tag/UMI locations and allowed mismatches.
3. Run the tag counter.

**Example config (`config.txt`):**
```txt
## Do not change the order of the configuration lines

# cell barcode location
-cb_location 0:17
-cb_mismatch 0

# sampletag barcode location
-sb_location 25:41
-sb_seq TCAGATGCCTGGATCC,AGCTGCTAACTGCAGC,ATTCAAGGGCAGCCGC,CACTAGGATCCGATGC
-sb_mismatch 1

# UMI location
-umi_location 17:29
```

**Run:**
```bash
python count_heatag_tags.py   -r1 file_R1.fq.gz   -r2 file_R2.fq.gz   -c config.txt   -b barcodes.tsv   -p mr   -t 10
```

> **Notes**
> - Positions use `start:end` with 0‑based indices; slicing semantics follow the script implementation.
> - `-sb_seq` is a comma‑separated list of HEATag barcodes.
> - Check `python count_heatag_tags.py -h` for all flags.

---

## 1) HEATag Barcode Counting

- **Script**: `count_heatag_tags.py`
- **Inputs**: `R1/R2 FASTQ`, `barcodes.tsv`, `config.txt`
- **Outputs**: per‑cell tag counts/assignments (TSV/CSV, per script defaults)

Configuration keys:
- `-cb_location`, `-cb_mismatch` – cell barcode slice and allowed mismatches
- `-sb_location`, `-sb_seq`, `-sb_mismatch` – sample‑tag slice, sequences, and allowed mismatches
- `-umi_location` – UMI slice

Tips:
- Keep the order of config lines unchanged (the parser expects it).
- Start with strict mismatch settings; relax only if necessary after QC.

---

## 2) Genome Functional Annotation

- `blast_uniprot.sh` — BLAST against local **UniProt** db  
- `add_uniprot_annotations.py` — add UniProt annotation to BLAST hits  
- `blast_nr.sh` — BLAST against local **NR** db  
- `reciprocal_best_hits.sh` — **Reciprocal Best Hits** pipeline to infer putative orthologs

> Ensure local BLAST databases are indexed and paths in the scripts are correct for your environment.

---

## 3) Single‑Cell RNA‑seq Processing

- `functions.R` — a collection of utility functions used across analyses  
- `scRNA_pipeline.R` — main pipeline for downstream scRNA‑seq analysis  
- `filter_cross_species_cells.R` — routines for removing cross‑species contaminants  
- `sample_demux.R` — sample demultiplexing used in this study  
- `differential_expression_GO.R` — differential expression and GO enrichment analysis

> Run these scripts in R (≥4.x). Each script contains its primary entry points.

---

## 4) SNP‑based Demultiplexing

- `scsplit_prep.sh` — prepare BAM indices, etc.  
- `scsplit_call_snvs.sh` — call single‑nucleotide variants (SNVs)  
- `freemuxlet_demux.sh` — demultiplex with **Freemuxlet**

---

## 5) Visualization

- `Visualization_funtion.R` — custom functions for drawing graphics  
- `Visualization.R` — visualization code  

## Cite

If you use this code, please cite the HEATag paper:

> *A universal cost‑efficient sample labeling approach for multiplexed single‑cell RNA‑seq based on recombinant HUH‑endonuclease–agglutinin tagging.*  
> (DOI once available)

---

## License
This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.
---

## Contact
For questions or issues, please open a GitHub issue in this repository.
