#!/bin/bash

DIAG_TSV=##Put filepath here
UNDIAG_TSV=##Put filepath here
OUTPUT=##Put filepath here
FOLDX_DIR=##Put directory of CSV files with FoldX scores with uniprot ID as filename

# Read through the file line by line, skipping the first line (header)
{
    # Assign the header to a variable and add gnomad columns, skip header in loop
    read -r header
    long_header="$header$(echo -e '\tddG_fold\tfold_rank\tabs_fold_rank')"
    echo -e "$long_header" > $OUTPUT
    # Process each row
    while IFS=$'\t' read -r gnomad_id gene uniprot_id codon uniprot_start wild_type mutant rest; do
        if [[ "$mutant" =~ ^[ACDEFGHIKLMNPQRSTVWY]$ && "$mutant" != "$wild_type" ]]; then
            foldx_line=$(grep "${wild_type}${uniprot_start}${mutant}" "${FOLDX_DIR}${uniprot_id}.csv")
            IFS=',' read -r variant ddG_fold fold_rank abs_fold_rank <<< "$foldx_line"
            if [[ "$ddG_fold" == "" ]]; then
                ddG_fold="."
            fi
       	    if [[ "$fold_rank" == "" ]]; then
       	       	fold_rank="."
       	    fi
       	    if [[ "$abs_fold_rank" == "" ]]; then
       	       	abs_fold_rank="."
       	    fi
        else
            ddG_fold="."
            fold_rank="."
            abs_fold_rank="."
        fi
        echo -e "$gnomad_id\t$gene\t$uniprot_id\t$codon\t$uniprot_start\t$wild_type\t$mutant\t$rest\t$ddG_fold\t$fold_rank\t$abs_fold_rank" >> $OUTPUT
    done
} < "$UNDIAG_TSV"
