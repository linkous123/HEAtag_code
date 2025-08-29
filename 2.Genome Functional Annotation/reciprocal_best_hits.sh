#!/bin/bash

for f in Amer_marker.fa Pver.protein.faa; do
  awk '/^>/{split($0,a," "); print a[1]; next} {print}' "$f" > "${f%.fa*}.clean.fa"
done

seqkit seq -u -w 0 Pver.protein.clean.fa \
| seqkit replace -s -p '[^ACDEFGHIKLMNPQRSTVWYBXZOU]' -r X \
> Pver.protein.diamond.fa

diamond makedb -p 8 --in Pver.protein.clean.fa -d Pver_db


diamond blastp -p 8 \
  -q Amer_marker.clean.fa -d Pver_db \
  -e 1e-5 --max-target-seqs 25 --very-sensitive \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen \
  -o amer2pver.m6

  awk 'BEGIN{FS=OFS="\t"}
     {
       qid=$1; sid=$2; pident=$3; alen=$4; e=$5; bit=$6; qlen=$7; slen=$8;
       qcov=alen/qlen; scov=alen/slen;
       pass=(pident>=30 && qcov>=0.50 && scov>=0.50 && e<=1e-5);
       if(pass && ( !(qid in best) || bit>bestbit[qid] )){best[qid]=$0; bestbit[qid]=bit}
     }
     END{for(q in best) print best[q]}' amer2pver.m6 \
| sort > amer2pver.best.tsv


cut -f2 amer2pver.best.tsv | sort -u > pver_hits.list
seqkit grep -n -f pver_hits.list Pver.protein.clean.fa > Pver_hits.fa

diamond makedb -p 8 --in Amer.protein.clean.fa -d Amer_db

diamond blastp -p 8 \
  -q Pver_hits.fa -d Amer_db \
  -e 1e-5 --max-target-seqs 1 --very-sensitive \
  --outfmt 6 qseqid sseqid bitscore \
  -o pver2amer.top1.m6

awk 'NR==FNR{rev[$1]=$2; next} {if(rev[$2]==$1) print $0"\tRBH"}' \
    pver2amer.top1.m6 amer2pver.best.tsv > amer_pver.RBH.tsv