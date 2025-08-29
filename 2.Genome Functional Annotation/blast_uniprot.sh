#!/bin/bash

blastp -query {input} -db db/SwissProt/uniprot_sprot -evalue 1e-5 -outfmt 6 -out {output} -max_target_seqs 6 -num_threads 32

awk '
        BEGIN {{ FS = "\t"; OFS = "\t" }}
        {{
            query = $1;
            if (query != last_query) {{
                if (last_query != "") {{
                    print best_result;
                }}
                last_query = query;
                max_identity = $3;
                best_result = $0;
            }} else {{
                if ($3 > max_identity) {{
                    max_identity = $3;
                    best_result = $0;
                }}
            }}
        }}
        END {{ if (last_query != "") print best_result }}
' {input} > {output}


python add_uniprot_annotations.py -b {input} -u db/SwissProt/fasta/uniprot_annotations.txt -o {output}