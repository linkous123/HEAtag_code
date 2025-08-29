#!/bin/bash

VCF="st_3testis/freebayes_var.vcf"
SCSPLIT_OUTDIR="st_3testis"
BARCODES="mtx/barcodes.tsv"
N="3"
FREEMUXLET_OUTDIR="freemuxlet"
UMI_TAG="MA"

singularity exec Demuxafy.sif scSplit count -c $VCF -v $SCSPLIT_OUTDIR/freebayes_var_qual30.recode.vcf -i $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam -b $BARCODES -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -o $SCSPLIT_OUTDIR

singularity exec Demuxafy.sif scSplit run -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -n $N -o $SCSPLIT_OUTDIR
# -n 0: autodetect mode

singularity exec Demuxafy.sif scSplit genotype -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -p $SCSPLIT_OUTDIR/scSplit_P_s_c.csv -o $SCSPLIT_OUTDIR


singularity exec Demuxafy.sif popscle_pileup.py \
--sam $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam \
--vcf $VCF \
--tag-UMI $UMI_TAG \
--group-list $BARCODES \
--out $FREEMUXLET_OUTDIR/pileup

singularity exec Demuxafy.sif popscle freemuxlet \
        --plp $FREEMUXLET_OUTDIR/pileup \
        --out $FREEMUXLET_OUTDIR/freemuxlet \
        --group-list $BARCODES \
        --nsample $N

singularity exec Demuxafy.sif bash Freemuxlet_summary.sh $FREEMUXLET_OUTDIR/freemuxlet.clust1.samples.gz

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
        -r $VCF \
        -c $FREEMUXLET_OUTDIR/freemuxlet.clust1.vcf.gz \
        -o $FREEMUXLET_OUTDIR