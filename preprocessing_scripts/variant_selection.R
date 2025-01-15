#!/bin/R
library(tidyverse)

# parse file names passed by the command line
args <- commandArgs(trailingOnly = TRUE)
vep_output <- args[1]

# fetch UniProt release version
uniprot_release <- read_lines(
  file = paste0(
    'https://ftp.uniprot.org/pub/databases/',
    'uniprot/current_release/relnotes.txt'
  ),
  n_max = 1
) |>
  str_remove('UniProt Release ')

# canonical and full UniProt proteomes
canonical_proteome <- paste0('UP000005640_', uniprot_release, '.fa.gz')
if (!file.exists(canonical_proteome)) {
  download.file(
    url = paste0(
      'https://ftp.uniprot.org/pub/databases/uniprot/',
      'current_release/knowledgebase/reference_proteomes/',
      'Eukaryota/UP000005640/UP000005640_9606.fasta.gz'
    ),
    destfile = canonical_proteome
  )
}

complete_proteome <- paste0('complete_', canonical_proteome)
if (!file.exists(complete_proteome)) {
  download.file(
    url = paste0(
      'https://rest.uniprot.org/uniprotkb/stream?format=fasta',
      '&query=%28%28proteome%3AUP000005640%29%29'
    ),
    destfile = complete_proteome
  )
}

# parse proteomes
canonical_fa <- read_lines(canonical_proteome)
complete_fa <- read_lines(complete_proteome)

pull_fa <- function(raw_fa) {
  # Parses a UniProt fasta as a named list, where names are accession IDs.
  #
  # Arguments:
  #   url (character): URL of the fasta file
  #
  # Returns:
  #   Named list whose elements are single character vectors containing
  #   the protein sequence. Names are UniProt accession IDs.

  header <- str_subset(raw_fa, '>')
  uniprot_acc <- str_extract(header, '(?<=\\|).+(?=\\|)')
  sep <- str_detect(raw_fa, '>')
  fa <- split(raw_fa[!sep], cumsum(c(TRUE, diff(sep) < 0))[!sep])
  fa <- map(fa, paste0, collapse = '')
  names(fa) <- uniprot_acc
  return(fa)
}

# store sequences as a named list
fa <- pull_fa(canonical_fa)
all_fa <- pull_fa(complete_fa)

# gene to canonical UniProt accession ID mapping table
gene_to_canonical_uniprot <- canonical_fa |>
  str_subset('>') |>
  as_tibble_col(column_name = 'header') |>
  transmute(gene = str_extract(header, '(?<=GN=)(.+)(?=\\ PE=)'),
            uniprot_id = str_extract(header, '(?<=\\|).+(?=\\|)'))

# extract all UniProt IDs and sequences with identical sequences to canonical
identical_seq_to_canonical <- complete_fa |>
  str_subset('>') |>
  as_tibble_col(column_name = 'header') |>
  transmute(gene = str_extract(header, '(?<=GN=)(.+)(?=\\ PE=)'),
            uniprot_id = str_extract(header, '(?<=\\|).+(?=\\|)')) |>
  mutate(seq = as.character(all_fa[uniprot_id])) |>
  add_count(gene, seq) |>
  filter(n > 1) |>
  filter(!uniprot_id %in% names(fa)) |>
  select(gene, secondary_id = uniprot_id, seq) |>
  inner_join(gene_to_canonical_uniprot, by = join_by(gene)) |>
  rename(primary_id = uniprot_id)

# append the canonical gene-UniProt mapping table
gene_to_canonical_uniprot <- bind_rows(
  gene_to_canonical_uniprot,
  select(identical_seq_to_canonical, gene, uniprot_id = secondary_id))

# append the canonical sequence list
fa <- append(fa, with(identical_seq_to_canonical, split(seq, secondary_id)))

# parse Ensembl VEP output
vep_mapping <- read_tsv(
  file = vep_output,
  skip = 44,
  col_select = c(
    id = `#Uploaded_variation`,
    codon = Codons,
    swissprot = SWISSPROT,
    trembl = TREMBL,
    isoform = UNIPROT_ISOFORM,
    uniprot_start = Protein_position,
    aa_change = Amino_acids,
    consequence = Consequence
  ),
  col_types = 'cccccccc'
) |>
  # exclude those without UniProt start position
  filter(uniprot_start != '-') |>
  mutate(
    # merge SwissProt IDs with TrEMBL IDs
    uniprot_id = if_else(swissprot == '-', trembl, swissprot),
    # remove version number if present
    uniprot_id = str_extract(uniprot_id, '[^.]+'),
    # extract isoform number
    isoform = str_extract(isoform, '(?<=-).*'),
    # parse wild-type and mutant amino acids
    wild_type = str_extract(aa_change, '[^/]+'),
    mutant = str_extract(aa_change, '(?<=/)[^/].*'),
    # for synonymous changes, leave wild type as mutant
    mutant = if_else(is.na(mutant), wild_type, mutant),
    # create UniProt end for indels
    uniprot_start = as.integer(parse_number(uniprot_start)),
    uniprot_end = uniprot_start + str_length(wild_type) - 1L
  ) |>
  # check residue correspondence to wild-type amino acid
  filter(wild_type == '*' |
           wild_type == str_sub(as.character(fa[uniprot_id]),
                                uniprot_start, uniprot_end)
  ) |>
  # select the isoform with the most number of variants
  add_count(uniprot_id, isoform) |>
  arrange(desc(n), str_detect(consequence, 'missense')) |>
  # keep distinct allele IDs
  distinct(id, .keep_all = TRUE) |>
  # add gene names
  inner_join(gene_to_canonical_uniprot, by = join_by(uniprot_id)) |>
  # convert secondary IDs with identical sequence to primary IDs
  left_join(select(identical_seq_to_canonical, primary_id, secondary_id),
            by = join_by(uniprot_id == secondary_id)) |>
  mutate(uniprot_id = if_else(is.na(primary_id), uniprot_id, primary_id)) |>
  select(id, gene, uniprot_id, codon, uniprot_start, wild_type, mutant, consequence)

# write out
write_tsv(x = vep_mapping,
          file = paste0(
            str_remove(vep_output, '_minimal_uniprot_raw.tsv'),
            '_uniprot_',
            uniprot_release,
            '.tsv.gz'
          ))

# done
