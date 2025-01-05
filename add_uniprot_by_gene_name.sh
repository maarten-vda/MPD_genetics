#!/bin/bash

# Check if the input file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file.tsv>"
    exit 1
fi

# Input file
input_file="$1"

# Check if the file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# Read the file line by line
while IFS=$'\t' read -r col1 rest; do
    col1=$(echo "$col1" | sed 's/ //g' | tr -d '\n')
    rest=$(echo "$rest" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')  # Trim leading and trailing spaces
    response=$(curl -sH "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606+AND+gene:${col1}")
    uniprot=$(echo "$response" | awk -F'\t' 'NR==2 {print $1}')
    if [[ $uniprot == "" ]]; then
        uniprot="."
    fi
    #echo $rest
    echo -e "$col1\t$rest\t$uniprot"
done < "$input_file"

