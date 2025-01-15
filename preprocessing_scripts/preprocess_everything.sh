
#!/bin/bash

#Assumes R, Ensembl VEP, and bcftools are in filepath


INPUT_VCF=## Put filepath here
MINIMAL_VCF_FILE=## Put filepath here
VEP_OUTPUT_FILE="variants_vep_output.tsv"

bcftools annotate -x INFO $INPUT_VCF | bcftools view -f PASS | bcftools norm -m -any -Oz -o $MINIMAL_VCF_FILE

fast_add_gnomad.sh $MINIMAL_VCF_FILE gnomad_annotated.tsv

# run Ensembl VEP
vep -i gnomad_annotated.tsv \
    -o "$VEP_OUTPUT_FILE" \
    --format "vcf" \
    --assembly GRCh38 \
    --cache --offline \
    --dir /exports/igmm/eddie/marsh-lab/data/ensembl_vep112/ \
    --coding_only --uniprot --no_stats --tab \
    --buffer_size 20000 --fork 4

Rscript variant_selection.R "$VEP_OUTPUT_FILE"

##Output name changes uniprot release

UNIPROT_RELEASE=##Uniprot release here

CODING_VARIANTS=$(echo "$VEP_OUTPUT_FILE" | sed 's/_minimal_uniprot_raw\.tsv$//')"_uniprot_${UNIPROT_RELEASE}.tsv.gz")

bash add_summary_stats.sh "$CODING_VARIANTS" coding_variants_variant_summaries.tsv
