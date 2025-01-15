## Assumes bcftools is in $PATH


# Input string
INPUT="$1"
gnomad_ID="$2"
PED"$3"
VCF="$4"


IFS='-' read -r CHR POS REF ALT <<< "$gnomad_ID"

# Find rows matching column 2
MATCHES=$(awk -v input="$INPUT" '$2 == input' "$PED")

read family proband father mother <<< $(echo -e "$MATCHES" | awk -F '\t' '{print $1, $2, $3, $4}')
probands=$(bcftools query -r "chr$CHR:$POS" -i "REF=\"$REF\" & ALT=\"$ALT\"" $VCF -f '[%SAMPLE:%GT\n]' -s "$proband,$father,$mother" --force-samples)
echo "#Family $family genotypes at position $gnomad_ID"
echo $probands
