#!/bin/bash

# Input and output file paths
INPUT=$1
OUTPUT=$2

# Process the file
# The command retains the header, then filters rows where column 55 is <= 0.0001
awk 'BEGIN { FS="\t"; OFS="\t" } \
    NR == 1 || ($55 <= 0.0001 && $55 != ".") { print }' "$INPUT" > "$OUTPUT"

# Completion message
echo "Filtered rows written to '$OUTPUT'."

