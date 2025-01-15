import argparse
import pandas as pd
import numpy as np
from queue import Queue

def get_inheritance(omim_tsv, uniprot_id):
    # Gets inheritance modes from OMIM TSV file
    try:
        inheritance = omim_tsv.loc[omim_tsv['uniprot_id'] == uniprot_id]['inheritance'].values
    except IndexError:
        inheritance = None
    return inheritance

def get_vep_scores(vep_scores, substitution):
    # Gets VEP scores from VEP directory
    try:
        correct_row = vep_scores.loc[vep_scores['variant'] == substitution]
        if "ESM-1v" in correct_row.columns:
            esm_1v = correct_row['ESM-1v'].values[0]
        else:
            esm_1v = None
        if "AlphaMissense" in correct_row.columns:
            alphamissense = correct_row['AlphaMissense'].values[0]
        else:
            alphamissense = None
        if "CPT" in correct_row.columns:
            cpt = correct_row['CPT'].values[0]
        else:
            cpt = None
        if "GEMME" in correct_row.columns:
            gemme = correct_row['GEMME'].values[0]
        else:
            gemme = None
        if "popEVE" in correct_row.columns:
            popeve = correct_row['popEVE'].values[0]
        else:
            popeve = None

    except Exception as e:
        print(e)
        esm_1v, alphamissense, cpt, gemme, popeve = None, None, None, None, None
    
    return(esm_1v, alphamissense, cpt, gemme, popeve)


def worker(queue, modified_rows, omim_tsv, vep_dir):
    while not queue.empty():
        index, row = queue.get()
        try:
            uniprot_id = row['uniprot_id']
            wild_type = row['wild_type']
            mutant = row['mutant']
            uniprot_start = row['uniprot_start']
            substitution = str(wild_type) + str(uniprot_start) + str(mutant)
            consequence = row['consequence']

            # Get inheritance information from the OMIM data
            inheritance = get_inheritance(omim_tsv, uniprot_id)

            # Initialize variables for VEP scores
            esm_1v, alphamissense, cpt, gemme, popeve = None, None, None, None, None
            
            # If the consequence is missense_variant, get VEP scores
            if consequence == "missense_variant":
                try:
                    vep_scores = pd.read_csv(f"{vep_dir}/{uniprot_id}.csv")
                    esm_1v, alphamissense, cpt, gemme, popeve = get_vep_scores(vep_scores, substitution)
                    print(esm_1v, alphamissense, cpt, gemme, popeve)
                except Exception as e: ### This line used to be filenotfounderror but I expanded it to all exceptions since I think theres some CSVs that are empty but exist
                    print(e)
                    pass  # Keep None values as initialized

            # Append the modified row with new columns
            modified_row = row.tolist() + [inheritance.tolist(), esm_1v, alphamissense, cpt, gemme, popeve]
            modified_rows.append(modified_row)

        finally:
            queue.task_done()

def main(probands, omim, vep_dir):
    # Read the input files
    proband_tsv = pd.read_csv(probands, sep='\t')
    omim_tsv = pd.read_csv(omim, sep='\t')
    
    # Define the header for the output DataFrame
    header = proband_tsv.columns.tolist() + ['inheritance', 'ESM-1v', 'AlphaMissense', 'CPT', 'GEMME', 'popEVE']
    
    # Create a list to store modified rows
    modified_rows = []
    
    # Iterate over each row in the probands DataFrame
    for index, row in proband_tsv.iterrows():
        uniprot_id = row['uniprot_id']
        wild_type = row['wild_type']
        mutant = row['mutant']
        uniprot_start = row['uniprot_start']
        substitution = str(wild_type) + str(uniprot_start) + str(mutant)
        consequence = row['consequence']
        
        # Get inheritance information from the OMIM data
        inheritance = get_inheritance(omim_tsv, uniprot_id)
        
        # Initialize variables for VEP scores
        esm_1v, alphamissense, cpt, gemme, popeve = None, None, None, None, None
        
        # If the consequence is missense_variant, get VEP scores
        if consequence == "missense_variant":
            try:
                vep_scores = pd.read_csv(f"{vep_dir}/{uniprot_id}.csv")
                esm_1v, alphamissense, cpt, gemme, popeve = get_vep_scores(vep_scores, substitution)
                print(esm_1v, alphamissense, cpt, gemme, popeve)
            except FileNotFoundError:
                print(f"File not found: {vep_dir}/{uniprot_id}.csv")
                # Keep None values as initialized
            except pd.errors.EmptyDataError:
                print(f"Empty data: {vep_dir}/{uniprot_id}.csv")
                # Keep None values as initialized
            except Exception as e:
                print(e)
                # Keep None values as initialized

        previous_uniprot_id = uniprot_id

        if esm_1v is None or esm_1v == np.nan:
            esm_1v = "."
        if alphamissense is None or alphamissense == np.nan:
            alphamissense = "."
        if cpt is None or cpt == np.nan:
            cpt = "."
        if gemme is None or gemme == np.nan:
            gemme = "."
        if popeve is None or popeve == np.nan:
            popeve = "."

        # Append the modified row with new columns
        modified_row = row.tolist() + [inheritance.tolist(), esm_1v, alphamissense, cpt, gemme, popeve]
        modified_rows.append(modified_row)
    
    # Create a new DataFrame with the modified rows and the updated header
    modified_df = pd.DataFrame(modified_rows, columns=header)
    
    # Output the DataFrame
    return modified_df 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prints the contents of a TSV file line by line.")
    parser.add_argument('proband_tsv', type=str, help="Path to the TSV file to be read")
    parser.add_argument('omim_tsv', type=str, help="Path to the TSV file to be read")
    parser.add_argument('vep_dir', type=str, help="Path to the directory of VEP scores")
    parser.add_argument('output', type=str, help="Path to the output file")

    args = parser.parse_args()
    main(args.proband_tsv, args.omim_tsv, args.vep_dir).to_csv(args.output, sep='\t', index=False)



