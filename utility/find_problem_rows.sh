#!/bin/bash

# Input file
input_file="$1"

# Get the number of columns in the header row
header=$(head -n 1 "$input_file")
header_columns=$(awk -F'\t' '{print NF}' <<< "$header")

# Initialize row counter
row_number=0

# Process the file line by line
while IFS=$'\n' read -r line; do
    # Increment row counter
    ((row_number++))

    # Count the number of columns in the current row
    num_columns=$(awk -F'\t' '{print NF}' <<< "$line")

    # Check if the number of columns matches the header
    if (( num_columns != header_columns )); then
        echo "Row $row_number does not have $header_columns columns:"
        echo "$line"
    fi
done < "$input_file"

