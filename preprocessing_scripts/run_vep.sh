

if [ -z "$1" ]; then
  echo "Usage: $0 <vcf_file>"
  exit 1
fi

VCF_FILE=$1
ROOT_NAME=$(basename "$VCF_FILE" | awk -F '.' '{print $1}')
VEP_OUTPUT_FILE="${ROOT_NAME}_uniprot_raw.tsv"

##Assumes bcftools and Ensembl VEP are in $PATH

# run Ensembl VEP
vep -i "$VCF_FILE" \
    -o "$VEP_OUTPUT_FILE" \
    --format "vcf" \
    --assembly GRCh38 \
    --cache --offline \
    --dir /exports/igmm/eddie/marsh-lab/data/ensembl_vep112/ \
    --coding_only --uniprot --no_stats --tab \
    --buffer_size 20000 --fork 4

exit 0
# done
