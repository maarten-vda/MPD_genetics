#SCRIPTS DIRECTORY FOR ANALYSIS SCRIPTS

-add_inheritance_and_vep_scores.py

   -->Adds inheritance modes from OMIM and VEP scores from VEP Directory, where each file is named a uniprot ID and has comma separated VEP predictions for each missense variant in uniprot.
   
-add_foldx_scores.sh

   -->Adds FoldX scores to cohort coding variant TSV. 
   
-format_proband_ids.py

   -->Adds which probands are de novo, inherited from one parent, both parents, or unknown.
   
-add_gnomad_mafs.sh/format_gnomad_info.py

   -->Retrieves variant information from gnomAD genomes file, and quickly adds it to TSV
   
-add_inheritance_info.sh/format_inheritance_info.py

   -->Quickly retrieves parental genotypes from VCF file, formats and adds them to TSV
   
--add_go_enrichment.py

   -->Adds per-gene per-ontology GO enrichment scores for the disease, requires a list of all genes annotated with enrichment scores, see go_enrichment_analysis folder
