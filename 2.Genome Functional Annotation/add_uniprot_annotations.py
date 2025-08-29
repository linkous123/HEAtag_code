import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Merge BLAST results with UniProt annotations.")
    parser.add_argument('-b', '--blast', required=True, help="Path to the BLAST result file.")
    parser.add_argument('-u', '--uniprot', required=True, help="Path to the UniProt annotations file.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output file.")
    return parser.parse_args()

def main():
    args = parse_args()

    blast_result = pd.read_csv(args.blast, sep='\t', header=None, names=[
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ])

    uniprot_annotations = pd.read_csv(args.uniprot, sep='\t')

    blast_result['sseqid'] = blast_result['sseqid'].astype(str)
    uniprot_annotations['UniProt_ID'] = uniprot_annotations['UniProt_ID'].astype(str)

    merged_result = pd.merge(blast_result, uniprot_annotations, left_on='sseqid', right_on='UniProt_ID', how='left')

    final_result = merged_result[[
        'qseqid', 'pident', 'evalue', 'bitscore', 'UniProt_ID', 'UniProt_acc', 
        'gene_name', 'protein', 'Organism', 'Organism_Taxonomy', 'PE', 'SV'
    ]]

    final_result.to_csv(args.output, sep=',', index=False, header=True)
    print(f"Output saved to {args.output}")

if __name__ == "__main__":
    main()

