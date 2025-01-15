# Scripts Directory for Preprocessing Scripts

## `preprocess_everything.sh`
- **Description:**
  - From the cohort VCF, this script:
    - Removes the info field.
    - Adds gnomAD-style IDs.
    - Filters for PASS variants.
    - Splits multiallelic sites.
    - Runs Ensembl VEP.
    - Executes the `variant_selection.R` script to add UniProt-level identifiers.
    - Adds cohort summary stats per variant.
  - **Purpose:** This script performs all the preprocessing tasks.

---

## `run_vep.sh`
- **Description:**
  - Runs Ensembl VEP.
  - Executes the `variant_selection.R` script to correctly format the data.

---

## `variant_selection.R`
- **Description:**
  - Adds UniProt-level identifiers to the dataset.

---

## `add_summary_stats_to_unprocessed_tsv.sh`
- **Description:**
  - Adds summary statistics to a cohort TSV, including:
    - Number of variants.
    - Homozygous probands.
    - Allele numbers (ANs) for both diagnosed and undiagnosed cohorts.

---

## `fast_add_gnomad.sh`
- **Description:**
  - Adds gnomAD-style IDs to the cohort VCF.
  - Ensures these IDs are preserved when Ensembl VEP converts the data to TSV format.
