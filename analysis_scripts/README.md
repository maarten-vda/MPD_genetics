# Scripts Directory for Analysis Scripts

## `add_inheritance_and_vep_scores.py`
- **Description:**
  - Adds inheritance modes from OMIM.
  - Retrieves VEP scores from the VEP directory, where:
    - Each file is named after a UniProt ID.
    - Files contain comma-separated VEP predictions for each missense variant in UniProt.

---

## `add_foldx_scores.sh`
- **Description:**
  - Adds FoldX scores to the cohort coding variant TSV.

---

## `format_proband_ids.py`
- **Description:**
  - Identifies and categorizes probands as:
    - De novo.
    - Inherited from one parent.
    - Inherited from both parents.
    - Unknown.

---

## `add_gnomad_mafs.sh` / `format_gnomad_info.py`
- **Description:**
  - Retrieves variant information from the gnomAD genomes file.
  - Quickly adds gnomAD data to the cohort TSV.

---

## `add_inheritance_info.sh` / `format_inheritance_info.py`
- **Description:**
  - Quickly retrieves parental genotypes from the VCF file.
  - Formats and adds the inheritance information to the TSV.

---

## `add_go_enrichment.py`
- **Description:**
  - Adds per-gene, per-ontology GO enrichment scores for the disease.
  - **Requirement:** A list of all genes annotated with enrichment scores (see the `go_enrichment_analysis` folder for details).
