#!/bin/bash

DIAG_TSV=## Put filepath here
UNDIAG_TSV=## Put filepath here
SCRIPT=## Put filepath here
DIAG_VCF=## Put filepath here
UNDIAG_VCF=## Put filepath here
OUTPUT=## Put filepath here
PED=## Put filepath here


# Read through the file line by line, skipping the first line (header)
{
    # Assign the header to a variable and add gnomad columns, skip header in loop
    read -r header
    long_header="$header$(echo -e '\tprobands\tn_monoallelic\tmonoallelic_probands\tn_biallelic\tbiallelic_probands\tn_denovo\tdenovo_probands\tn_unknown_inheritance\tunknown_inheritance_probands')"
    echo -e "$long_header" > $OUTPUT
    # Process each row
    while IFS=$'\t' read -r gnomad_id rest; do
        IFS="-" read -r CHR POS REF ALT <<< "$gnomad_id"
        probands=$(bcftools query -r "chr$CHR:$POS" -i "REF=\"$REF\" & ALT=\"$ALT\"" $UNDIAG_VCF -f '[%SAMPLE:%GT\n]')
        python $SCRIPT "$gnomad_id" "$rest" "$probands" --ped $PED >> $OUTPUT
    done
} < "$UNDIAG_TSV"
