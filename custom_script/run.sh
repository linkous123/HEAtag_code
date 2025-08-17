#!/bin/bash

#SBATCH -J rebuild_4
#SBATCH -n 50
#SBATCH --mem 100G
#SBATCH -p nol
#SBATCH -o rebuild_4.%j.out
#SBATCH -e rebuild_4.%j.err

python rebuild_3.py -r1 /Share/user/limaor/rawdata/1sampletag/20250721/testis-T4-unlable-STL3_L1_1.fq.gz -r2 /Share/user/limaor/rawdata/1sampletag/20250721/testis-T4-unlable-STL3_L1_2.fq.gz -c config.txt -b barcodes.tsv -p mr -t 50


scontrol show job $SLURM_JOB_ID
