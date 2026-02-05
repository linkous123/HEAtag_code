#!/bin/bash

BAM="./bam/merged.bam"
SCSPLIT_OUTDIR="./outs"

singularity exec Demuxafy.sif samtools view -@ 32 -b -S -q 10 -F 3844 $BAM > $SCSPLIT_OUTDIR/filtered_bam.bam

singularity exec Demuxafy.sif samtools rmdup $SCSPLIT_OUTDIR/filtered_bam.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam

singularity exec Demuxafy.sif samtools sort -@ 32 -o $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam

singularity exec Demuxafy.sif samtools index $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam
