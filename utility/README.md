# Utility Scripts for TSV and VCF Processing

This directory contains a collection of utility scripts designed to assist with processing TSV and VCF files. Below is a description of each script and its functionality.

---

## `add_uniprot_by_gene_name.sh`
- **Description:** 
  - Adds a column with UniProt IDs to a TSV file.
  - Assumes the first column of the input TSV contains gene names.
- **Usage:**
  - Ensure the TSV file has gene names in the first column before running this script.

---

## `find_problem_rows.sh`
- **Description:** 
  - Checks the input TSV file for rows with fewer columns than the header row.
  - Useful for identifying formatting issues or corrupted rows in large datasets.
- **Usage:**
  - Provide a TSV file as input to the script to perform the check.
  - ```bash bash add_uniprot_by_gene_name.sh <path_to_input_file.tsv>
---

## `get_parents_genotypes.sh`
- **Description:** 
  - Extracts parental genotypes from a VCF file.
  - Outputs the relevant genotypes for further analysis or downstream processing.
- **Usage:**
  - Input a valid VCF file containing parental genotype information.

---

## `hgnc_to_uniprot.py`
- **Description:** 
  - Converts HGNC identifiers to UniProt IDs.
  - Requires the `--input` variable to specify the HGNC file or identifier.
- **Usage:**
  - Pass an HGNC file or individual identifier to the script using the `--input` option.

---

### General Notes
- **Dependencies:** Ensure all required dependencies and permissions are met for each script. Some may require specific command-line tools or Python libraries.
- **File Compatibility:** These scripts are designed for TSV and VCF files. Ensure your files are properly formatted before using them.
- **Error Handling:** Use the `find_problem_rows.sh` script to verify your TSV files for formatting issues.

Feel free to reach out if you encounter issues or have questions regarding these scripts.
