import argparse
import sys
import pandas as pd
import numpy as np
import json
from scipy.stats import rankdata

'''
This script makes 3 TSVs with per-gene scores of GO enrichment to known MPD disease genes.
The score is equal to the sum of the weights of the GO terms associated with the gene normalised by number of GO terms.
Weights are equal to the count of GO terms in known MPD disease genes scaled by the number of probands where the gene is causative, normalised by the number of genes or information content. 

Inputs for this are the scaled counts for all GO terms in known MPD disease genes (F_counts_scaled_n_probands.tsv, P_counts_scaled_n_probands.tsv, C_counts_scaled_n_probands.tsv), 
goa_human.gaf (https://current.geneontology.org/annotations/goa_human.gaf.gz), and go.json (https://purl.obolibrary.org/obo/go.json)

Output is a TSV with columns uniprot_id, F_enrichment, P_enrichment, C_enrichment

'''


### Preloading functions
'''
These preloading functions speed up the processing ~1000x by preloading data into memory
'''
def preload_goa(goa_filepath):
    with open(goa_filepath, 'r') as file:
    # Filter out lines starting with '!'
        filtered_lines = [line for line in file if not line.startswith('!')]
            # Split filtered lines into columns
        data = [line.strip().split('\t') for line in filtered_lines]
        data = pd.DataFrame(data[1:])
        return (data)
    

def preprocess_json(json_file):
    """
    Preprocess the JSON file into a dictionary mapping IDs to their 'lbl' values.
    
    :param json_file: Path to the JSON file
    :return: A dictionary {id: lbl}
    """
    with open(json_file, 'r') as file:
        data = json.load(file)
    
    id_to_lbl = {}
    for graph in data.get("graphs", []):
        for node in graph.get("nodes", []):
            node_id = node.get("id")
            go_id = node_id.split("/")[-1]
            go_id = go_id.replace("_", ":")
            lbl = node.get("lbl")
            if go_id and lbl:
                id_to_lbl[go_id] = lbl

    return id_to_lbl


### Helper functions

def get_go_terms(goa_data, uniprot_id):

    '''
    This function gets the GO terms for a given uniprot ID.
    '''

    filtered_data = goa_data[goa_data.iloc[:, 1] == uniprot_id]

    go_f = filtered_data[filtered_data.iloc[:, 8] == "F"]
    go_p = filtered_data[filtered_data.iloc[:, 8] == "P"]
    go_c = filtered_data[filtered_data.iloc[:, 8] == "C"]

    go_f_id = list(set(go_f[4].tolist()))
    go_p_id = list(set(go_p[4].tolist()))
    go_c_id = list(set(go_c[4].tolist()))

    go_dict = {"F": go_f_id, "P": go_p_id, "C": go_c_id}    

    return (go_dict)



def translate_go_terms(goa_dict, go_dict):
    '''
    This function adds the GO term labels to the dictionary of all GO IDs.
    '''
    goa_dict["F_lbl"] = []
    goa_dict["P_lbl"] = []
    goa_dict["C_lbl"] = []
    for key, value in goa_dict.items():
        if key == "F":
            for go_id in value:
                go_lbl = go_dict[go_id]
                goa_dict["F_lbl"].append(go_lbl)
        elif key == "P":
            for go_id in value:
                go_lbl = go_dict[go_id]
                goa_dict["P_lbl"].append(go_lbl)
        elif key == "C":
            for go_id in value:
                go_lbl = go_dict[go_id]
                goa_dict["C_lbl"].append(go_lbl)
    return goa_dict

def uniprot_to_go_list(goa_data, go_json_data):
    '''
    This function creates a dictionary of uniprot IDs to GO terms.
    '''
    uniprot_to_go = {}
    for uniprot_id in goa_data.iloc[:, 1]:
        go_dict = get_go_terms(goa_data, uniprot_id)
        uniprot_to_go[uniprot_id] = go_dict
    return uniprot_to_go


def get_go_terms_optimized(goa_data):
    """
    Optimized function to create a dictionary of UniProt IDs to GO terms.
    """
    grouped = goa_data.groupby(goa_data.columns[1])  # Group by UniProt ID (column 1)
    uniprot_go_dict = {}
    for uniprot_id, group in grouped:
        go_f = group[group.iloc[:, 8] == "F"].iloc[:, 4].unique().tolist()
        go_p = group[group.iloc[:, 8] == "P"].iloc[:, 4].unique().tolist()
        go_c = group[group.iloc[:, 8] == "C"].iloc[:, 4].unique().tolist()
        uniprot_go_dict[uniprot_id] = {"F": go_f, "P": go_p, "C": go_c}
    return uniprot_go_dict


def translate_go_terms_optimized(uniprot_go_dict, go_json_data):
    """
    Optimized function to add GO term labels to the dictionary.
    """
    for uniprot_id, go_dict in uniprot_go_dict.items():
        for category in ["F", "P", "C"]:
            go_ids = go_dict[category]
            go_labels = [go_json_data.get(go_id, "Unknown") for go_id in go_ids]
            go_dict[f"{category}_lbl"] = go_labels
    return uniprot_go_dict


def uniprot_to_go_list_optimized(goa_data, go_json_data):
    """
    Create a dictionary of UniProt IDs to GO terms using optimized grouping.
    """
    uniprot_go_dict = get_go_terms_optimized(goa_data)
    uniprot_go_dict = translate_go_terms_optimized(uniprot_go_dict, go_json_data)
    return uniprot_go_dict


def get_go_scaled_count(score_df, go_id_list, global_freqs, go_set):

    '''
    This function takes a GO ID, gets the score from the score_df.
    '''

    cumulative_score = 0
    for go_id in go_id_list:
        row = score_df[score_df.iloc[:, 0] == go_id]
        if not row.empty:
            scaled_count_value = row['Count'].iloc[0]
            background_freq = global_freqs[go_set][go_id]
            weighted_score = scaled_count_value / background_freq
        elif row.empty:
            weighted_score = 0
        cumulative_score += weighted_score

    return cumulative_score


def get_global_counts(uniprot_go_dict):
    '''
    This function takes the uniprot_go_dict outputted by uniprot_to_go_list_optimized and returns the global frequency of GO terms split by F, P and C.
    '''

    global_counts = {"F": {}, "P": {}, "C": {}}
    uniprot_count = len(uniprot_go_dict)
    for uniprot_id, go_dict in uniprot_go_dict.items():
        for go_lbl, go_list in go_dict.items():
            if go_lbl in ["F", "P", "C"]:
                for go_id in go_list:
                    if go_id in global_counts[go_lbl]:
                        global_counts[go_lbl][go_id] += 1
                    else:
                        global_counts[go_lbl][go_id] = 1

    # Normalise the counts by the number of genes
    for go_lbl in global_counts.keys():
        for go_id in global_counts[go_lbl].keys():
            global_counts[go_lbl][go_id] = global_counts[go_lbl][go_id] / uniprot_count

    return global_counts


def flatten_dict_to_tsv(nested_dict, output_filepath):
    # Prepare a list to store the rows of the TSV
    rows = []
    
    # Iterate through each outer dictionary key and its associated nested dictionary
    for outer_key, inner_dict in nested_dict.items():
        row = [outer_key]  # Start the row with the outer key (first column)
        
        # Add values from the inner dictionary to the row
        for inner_key in inner_dict:
            row.append(inner_dict[inner_key])
        
        # Add the row to the list of rows
        rows.append(row)
    
    # Convert the rows into a DataFrame
    df = pd.DataFrame(rows, columns=['Key'] + list(next(iter(nested_dict.values())).keys()))
    

    # Calculate percentile ranks for F, P, and C enrichment scores
    df['F_enrichment_norm'] = rankdata(df['F_enrichment'], method='min') / len(df)  # Normalized F enrichment
    df['P_enrichment_norm'] = rankdata(df['P_enrichment'], method='min') / len(df)  # Normalized P enrichment
    df['C_enrichment_norm'] = rankdata(df['C_enrichment'], method='min') / len(df)  # Normalized C enrichment


    # Write the DataFrame to a TSV file
    df.to_csv(output_filepath, sep='\t', index=False)


## Main functions


def make_enrichment_scores_tsv(P_filepath, F_filepath, C_filepath, output_filepath, goa_filepath, go_json_filepath):
    # This needs to iterate through all genes in GO
    goa_data = preload_goa(goa_filepath)
    go_json_data = preprocess_json(go_json_filepath)
    P_data = pd.read_csv(P_filepath, sep='\t')
    F_data = pd.read_csv(F_filepath, sep='\t')
    C_data = pd.read_csv(C_filepath, sep='\t')

    uniprot_go_dict = uniprot_to_go_list_optimized(goa_data, go_json_data)
    global_freqs = get_global_counts(uniprot_go_dict)


    # Iterate through each uniprot_id in the dictionary and calculate the enrichment score for each of P, F, C

    total_uniprot_ids = len(uniprot_go_dict)

    for index, key in enumerate(uniprot_go_dict.keys()):  # Iterate through all uniprot IDs in GO
        print(f"Processing {index + 1} of {total_uniprot_ids} uniprot IDs.")
        for go_lbl, go_id_list in uniprot_go_dict[key].items(): # Iterate through all GO terms for the uniprot ID
            num_go_terms = len(go_id_list)
            if go_lbl == "F":
                score = get_go_scaled_count(F_data, go_id_list, global_freqs, "F")
                if num_go_terms > 0:
                    F_enrichment = score / num_go_terms
                else:
                    F_enrichment = 0    
            elif go_lbl == "P":
                score = get_go_scaled_count(P_data, go_id_list, global_freqs, "P")
                if num_go_terms > 0:
                    P_enrichment = score / num_go_terms
                else:
                    P_enrichment = 0
            elif go_lbl == "C":
                score = get_go_scaled_count(C_data, go_id_list, global_freqs, "C")
                if num_go_terms > 0:
                    C_enrichment = score / num_go_terms
                else:
                    C_enrichment = 0
        uniprot_go_dict[key]["F_enrichment"] = F_enrichment
        uniprot_go_dict[key]["P_enrichment"] = P_enrichment
        uniprot_go_dict[key]["C_enrichment"] = C_enrichment
    
    # Flatten the dictionary to a TSV
    flatten_dict_to_tsv(uniprot_go_dict, output_filepath)


    return(print("Saved enrichment scores to TSV."))
    
    

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process gnomAD data and convert a semicolon-delimited string into a dictionary.")
    parser.add_argument("--P", type=str, help="The GO Process annotation scaled counts for known MPD disease genes.")
    parser.add_argument("--F", type=str, help="The GO Function annotation scaled counts for known MPD disease genes.")
    parser.add_argument("--C", type=str, help="The GO Component annotation scaled counts for known MPD disease genes.")
    parser.add_argument("--output", type=str, help="Path to output TSV.")
    parser.add_argument("--goa", type=str, help="The path to the GO annotations file.")
    parser.add_argument("--go", type=str, help="Path to the GO terms JSON file.")

    
    # Parse the arguments
    args = parser.parse_args()
    make_enrichment_scores_tsv(args.P, args.F, args.C, args.output, args.goa, args.go)
    

if __name__ == "__main__":
    main()
