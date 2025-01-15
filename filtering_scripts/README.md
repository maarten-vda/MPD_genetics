# Scripts Directory for Filtering Scripts

## `AF_grpmax_filter.sh`
- **Description:**
  - Applies a hard filter based on `AF_grpmax`, which is the maximum allele frequency in any population in gnomAD.

---

## `filter_FS_MQ_QD.py`
- **Description:**
  - Filters variants based on the following metrics:
    - **FS (Fisher Strand Bias):** Measures strand bias in sequencing data.
    - **MQ (Mapping Quality):** Indicates the quality of read alignment.
    - **QD (Quality by Depth):** Represents QUAL/DP, a measure of quality normalized by depth.

---

## `top5_percent_vep_go_filter.py`
- **Description:**
  - Filters out variants that:
    - Are not in the top 5% of any VEP scores.
    - Are not in the top 25% of GO enrichment scores for each category. (~11.5% of all human genes have GO enrichment scores in top 25% of each category)
  - **Note:** This filter removes 88.5% of genes.
