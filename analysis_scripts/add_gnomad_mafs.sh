#!/bin/bash

### This script assumes bcftools and python 3.10 are in PATH

DIAG_TSV=## Put filepath here
UNDIAG_TSV= ## Put filepath here
SCRIPT=## Put filepath here
GNOMAD=## Put filepath here
OUTPUT=## Put filepath here

# Read through the file line by line, skipping the first line (header)
{
    # Assign the header to a variable and add gnomad columns, skip header in loop
    read -r header
    long_header="$header$(echo -e '\tAF\tAF_XX\tAF_XY\tAF_afr_XX\tAF_afr_XY\tAF_afr\tAF_ami_XX\tAF_ami_XY\tAF_ami\tAF_amr_XX\tAF_amr_XY\tAF_amr\tAF_asj_XX\tAF_asj_XY\tAF_asj\tAF_eas_XX\tAF_eas_XY\tAF_eas\tAF_fin_XX\tAF_fin_XY\tAF_fin\tAF_mid_XX\tAF_mid_XY\tAF_mid\tAF_nfe_XX\tAF_nfe_XY\tAF_nfe\tAF_raw\tAF_remaining_XX\tAF_remaining_XY\tAF_remaining\tAF_sas_XX\tAF_sas_XY\tAF_sas\tAF_grpmax\tFS\tMQ\tQD\tinbreeding_coeff\tspliceai_ds_max\tphylop')"
    echo -e "$long_header" > $OUTPUT
    # Process each row
    while IFS=$'\t' read -r gnomad_id rest; do
        IFS="-" read -r CHR POS REF ALT <<< "$gnomad_id"
        gnomad_query=$(bcftools query -r "chr$CHR:$POS" -i "REF=\"$REF\" && ALT=\"$ALT\"" -f '%INFO\n' $GNOMAD)
        python $SCRIPT "$gnomad_id" "$rest" "$gnomad_query" >> $OUTPUT
    done
} < "$DIAG_TSV"
