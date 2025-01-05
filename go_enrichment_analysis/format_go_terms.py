import argparse
import sys
import pandas as pd
import numpy as np
import json

'''
This script makes TSVs with the number of times a go term is found in known MPD disease genes 
scaled by the number of probands where the gene is causative.
Input for this is genes_found_protein_only.tsv, goa_human.gaf (https://current.geneontology.org/annotations/goa_human.gaf.gz), and go.json (https://purl.obolibrary.org/obo/go.json)
'''


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


def get_go_terms(goa_data, uniprot_id):

    filtered_data = goa_data[goa_data.iloc[:, 1] == uniprot_id]

    go_f = filtered_data[filtered_data.iloc[:, 8] == "F"]
    go_p = filtered_data[filtered_data.iloc[:, 8] == "P"]
    go_c = filtered_data[filtered_data.iloc[:, 8] == "C"]

    go_f_id = list(set(go_f[4].tolist()))
    go_p_id = list(set(go_p[4].tolist()))
    go_c_id = list(set(go_c[4].tolist()))

    go_dict = {"F": go_f_id, "P": go_p_id, "C": go_c_id}    

    return (go_dict)



def add_go_to_tsv(goa_filepath, tsv_file, go_json_filepath, output_file):
    # Iterate over the tsv file and add the GO terms to the last columns
    goa_data = preload_goa(goa_filepath)
    go_json_data = preprocess_json(go_json_filepath)
    data = pd.read_csv(tsv_file, sep='\t')
    
    # Iterate through each row in the DataFrame
    for index, row in data.iterrows():
        uniprot_id = row[2] 

        # Get the GO terms dictionary for the current uniprot_id
        goa_dict = get_go_terms(goa_data, uniprot_id)
        
        # Translate the GO terms into labels using the provided function
        translated_go_terms_dict = translate_go_terms(goa_dict, go_json_data)
        
        # Extract values for the keys "F", "P", "C", "F_lbl", "P_lbl", "C_lbl"
        data.at[index, 'F'] = str(translated_go_terms_dict["F"])
        data.at[index, 'P'] = str(translated_go_terms_dict["P"])
        data.at[index, 'C'] = str(translated_go_terms_dict["C"])
        data.at[index, 'F_lbl'] = str(translated_go_terms_dict["F_lbl"])
        data.at[index, 'P_lbl'] = str(translated_go_terms_dict["P_lbl"])
        data.at[index, 'C_lbl'] = str(translated_go_terms_dict["C_lbl"])

    # Save the updated DataFrame to a new file (or overwrite the original file)
    #data.to_csv(output_file, sep='\t', index=False) ## Uncomment this and next to save the GO terms file
    #print("GO terms written to output file.")

    F_counts = {}
    P_counts = {}
    C_counts = {}
    for index, row in data.iterrows():
        F = row['F']
        P = row['P']
        C = row['C']
        n_probands = row[1]
        F = F.replace("[", "").replace("]", "").replace("'", "").split(", ")
        P = P.replace("[", "").replace("]", "").replace("'", "").split(", ")
        C = C.replace("[", "").replace("]", "").replace("'", "").split(", ")
        ## This is when its weighted by the number of probands with the variant in that gene, make it +=1 for just per gene
        for go_id in F:
            if go_id in F_counts:
                F_counts[go_id] += n_probands
            else:
                F_counts[go_id] = n_probands
        for go_id in P:
            if go_id in P_counts:
                P_counts[go_id] += n_probands
            else:
                P_counts[go_id] = n_probands
        for go_id in C:
            if go_id in C_counts:
                C_counts[go_id] += n_probands
            else:
                C_counts[go_id] = n_probands
    

    ## Make the summed counts into a DataFrame
    F_counts_df = pd.DataFrame(list(F_counts.items()), columns = ['GO_ID', 'Count'])
    F_new_column_values = F_counts_df.iloc[:, 0].map(go_json_data)
    F_counts_df.insert(1, 'GO_Term', F_new_column_values)
    F_counts_df.to_csv("F_counts_weighted_n_probands.tsv", sep='\t', index=False)

    P_counts_df = pd.DataFrame(list(P_counts.items()), columns = ['GO_ID', 'Count'])
    P_new_column_values = P_counts_df.iloc[:, 0].map(go_json_data)
    P_counts_df.insert(1, 'GO_Term', P_new_column_values)
    P_counts_df.to_csv("P_counts_weighted_n_probands.tsv", sep='\t', index=False)

    C_counts_df = pd.DataFrame(list(C_counts.items()), columns = ['GO_ID', 'Count'])
    C_new_column_values = C_counts_df.iloc[:, 0].map(go_json_data)
    C_counts_df.insert(1, 'GO_Term', C_new_column_values)
    C_counts_df.to_csv("C_counts_weighted_n_probands.tsv", sep='\t', index=False)


    ## Next step - weight the counts by the number of probands with the variant in that gene

    return (print("Cumulative GO terms written to output files."))



def translate_go_terms(goa_dict, go_dict):
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


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process gnomAD data and convert a semicolon-delimited string into a dictionary.")
    parser.add_argument("--input", type=str, help="Path to genes_found_protein_only.tsv input file.")
    parser.add_argument("--output", type=str, help="Path to output.")
    parser.add_argument("--goa", type=str, help="The path to the GO annotations file.")
    parser.add_argument("--go", type=str, help="Path to the GO terms JSON file.")

    
    # Parse the arguments
    args = parser.parse_args()
    add_go_to_tsv(args.goa, args.input, args.go, args.output)
    

if __name__ == "__main__":
    main()
