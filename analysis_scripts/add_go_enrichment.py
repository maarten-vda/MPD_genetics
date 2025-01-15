import csv

# File paths
DIAG_TSV = ## Put filepath here
UNDIAG_TSV = ## Put filepath here
UNDIAG_OUT = ## Put filepath here
DIAG_OUT = ## Put filepath here
GO_TSV = ## Put filepath here

# Step 1: Load GO_TSV into a dictionary with uniprot_id as key
go_dict = {}
with open(GO_TSV, 'r') as go_file:
    go_reader = csv.reader(go_file, delimiter='\t')
    next(go_reader)  # Skip the header
    for row in go_reader:
        uniprot_id = row[0]
        F_score = row[10]  # Column 11
        P_score = row[11]  # Column 12
        C_score = row[12]  # Column 13
        go_dict[uniprot_id] = {'F': F_score, 'P': P_score, 'C': C_score}

# Step 2: Open the undiagnosed TSV file, process each line, and write output efficiently
with open(DIAG_TSV, 'r') as undiagnosed_file, open(DIAG_OUT, 'w', newline='') as undiagnosed_out:
    undiagnosed_reader = csv.reader(undiagnosed_file, delimiter='\t')
    undiagnosed_writer = csv.writer(undiagnosed_out, delimiter='\t')

    # Read header and add the new columns for F, P, C enrichment
    header = next(undiagnosed_reader)
    header.extend(['F_enrichment', 'P_enrichment', 'C_enrichment'])
    undiagnosed_writer.writerow(header)  # Write the updated header to output file

    # Step 3: Process each row from the undiagnosed file
    for row in undiagnosed_reader:
        gnomad_id, gene, uniprot_id, *rest = row

        # Retrieve GO scores from the dictionary (or 'NA' if not found)
        enrichment_scores = go_dict.get(uniprot_id, {'F': '.', 'P': '.', 'C': '.'})
        F_score = enrichment_scores['F']
        P_score = enrichment_scores['P']
        C_score = enrichment_scores['C']

        # Write the processed row with enrichment scores to the output file
        undiagnosed_writer.writerow([gnomad_id, gene, uniprot_id] + rest + [F_score, P_score, C_score])

