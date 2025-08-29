#!/bin/bash


diamond blastp -p 32 --query pv_clean.pep --db db/nr_diamond/nr.dmnd -e 1e-5 -k 1 --sensitive --outfmt 6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore --out pv_nr.result
