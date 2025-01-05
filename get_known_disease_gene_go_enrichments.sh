#!/bin/bash

# Path to the input files
genes_found_go_terms="genes_found_go_terms.tsv"
all_genes_go_enrichments="all_genes_go_enrichments.tsv"

# Iterate through each line in genes_found_go_terms (ignoring header)
tail -n +1 "$genes_found_go_terms" | while IFS=$'\t' read -r col1 col2 col3 rest; do
    # Search for the value of col3 (UniProt ID) in the first column of all_genes_go_enrichments
    # Print columns 11-13 of the matching row
    awk -F'\t' -v uniprot_id="$col3" '$1 == uniprot_id {print $11, $12, $13}' "$all_genes_go_enrichments" | while read var11 var12 var13; do
        echo -e "$col1\t$col2\t$col3\t$rest\t$var11\t$var12\t$var13"
    done
done

