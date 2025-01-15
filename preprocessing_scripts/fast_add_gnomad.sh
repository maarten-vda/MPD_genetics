
#!/bin/bash
##Assumes bcftools is in the $PATH

INPUT_FILE=$1
OUTPUT_FILE=$2

# Use awk to process the file in one go
awk '
  BEGIN { OFS="\t" }
  # Function to construct the gnomad_id
  function construct_gnomad_id(chr, pos, ref, alt) {
    gsub("chr", "", chr)
    gsub("\\*", "del", alt)
    return chr "-" pos "-" ref "-" alt
  }
  {
    # Skip lines starting with ##
    if ($0 ~ /^##/) {
      next
    }
    # Process lines starting with a single #
    else if ($0 ~ /^#/) {
      print $0 >> "'"${OUTPUT_FILE}"'"
    }
    # Process all other lines
    else {
      gnomad_id = construct_gnomad_id($1, $2, $4, $5)
      $3 = gnomad_id
      print $0 >> "'"${OUTPUT_FILE}"'"
    }
  }
' "$INPUT_FILE"

