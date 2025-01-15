# Utility Scripts for TSV and VCF Processing

This directory contains a collection of utility scripts designed to assist with processing TSV and VCF files. Below is a description of each script and its functionality.

---

## `add_uniprot_by_gene_name.sh`
- **Description:** 
  - Adds a column with UniProt IDs to a TSV file.
  - Assumes the first column of the input TSV contains gene names.
- **Usage:**
  - ```bash add_uniprot_by_gene_name.sh <path_to_input_file.tsv>```
  - Provide a TSV file as an argument where the first column contains gene names.

---

## `find_problem_rows.sh`
- **Description:** 
  - Checks the input TSV file for rows with fewer columns than the header row.
  - Useful for identifying formatting issues or corrupted rows in large datasets.
- **Usage:**
  - ```bash find_problem_rows.sh <path_to_input_file.tsv>```
  - Provide a TSV file as input to check for problem rows with mismatched columns.

---

## `get_parents_genotypes.sh`
- **Description:** 
  - Extracts parental genotypes from a VCF file.
  - Outputs the relevant genotypes for further analysis or downstream processing.
- **Usage:**
  - ```bash get_parents_genotypes.sh <VCF_file> <gnomad_ID> <PED_file> <output_file>```
  - The script requires the following inputs:
    - **VCF_file**: Path to the VCF file containing genotype data.
    - **gnomad_ID**: gnomAD identifier.
    - **PED_file**: Pedigree file (if applicable).
    - **output_file**: Path to store the output.

---

## `hgnc_to_uniprot.py`
- **Description:** 
  - Converts HGNC identifiers to UniProt IDs.
  - Requires the `--input` variable to specify the HGNC file or identifier.
- **Usage:**
  - ```python hgnc_to_uniprot.py --input <hgnc_file_or_identifier>```
  - Provide a file or individual identifier to convert from HGNC to UniProt IDs.

---

### General Notes
- **Dependencies:** Ensure all required dependencies and permissions are met for each script. Some may require specific command-line tools or Python libraries.
- **File Compatibility:** These scripts are designed for TSV and VCF files. Ensure your files are properly formatted before using them.
- **Error Handling:** Use the `find_problem_rows.sh` script to verify your TSV files for formatting issues.
