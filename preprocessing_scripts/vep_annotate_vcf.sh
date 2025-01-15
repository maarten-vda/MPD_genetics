#!/bin/sh

INPUT=## Put filepath here
OUTPUT=## Put filepath here

##Assumes bcftools and ensembl VEP are in $PATH

# run Ensembl VEP
vep -i "$INPUT" \
    -o "$OUTPUT" \
    --vcf \
    --assembly GRCh38 \
    --cache --offline \
    --dir /exports/igmm/eddie/marsh-lab/data/ensembl_vep112/ \
    --buffer_size 20000 --fork 4 --everything
