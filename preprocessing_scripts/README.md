# SCRIPTS DIRECTORY FOR PREPROCESSING SCRIPTS

--preprocess_everything.sh

  -->From the cohort VCF, this script removes the info field, adds gnomAD-style IDs, filters for pass variants, splits multiallelic sites, runs Ensembl VEP, and then the variant_selection.R script adding uniprot-level identifiers, and adds cohort summary stats per variant. 
  
  -->This script does all the preprocessing 

--run_vep.sh

  -->Runs ensembl VEP and then the variant selection R script to correctly format the data



--variant_selection.R

  -->Adds uniprot-level identifiers



--add_summary_stats_to_unprocessed_tsv.sh

  -->Adds summary statistics for how many variants, homozygous probands, and ANs for both diagnosed and undiagnosed cohorts to a cohort TSV.



--fast_add_gnomad.sh

  -->This adds gnomAD style IDs to the cohort VCF, these are preserved when ensembl VEP converts to TSV. 

