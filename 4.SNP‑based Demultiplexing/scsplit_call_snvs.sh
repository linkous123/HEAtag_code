#!/bin/bash

FASTA="/path/GRCm39.genome.fa"
SCSPLIT_OUTDIR="t34"


singularity exec Demuxafy.sif freebayes -f $FASTA -iXu -C 2 -q 1 $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam > $SCSPLIT_OUTDIR/freebayes_var.vcf

singularity exec Demuxafy.sif vcftools --gzvcf $SCSPLIT_OUTDIR/freebayes_var.vcf --minQ 30 --recode --recode-INFO-all --out $SCSPLIT_OUTDIR/freebayes_var_qual30
