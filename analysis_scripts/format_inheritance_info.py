import argparse
import sys
import pandas as pd

def probands_to_dict(probands):
    '''Iterate through the lines of the proband string, separate each line by the : and make the part before : the key and the part after : the value'''
    probands_dict = {}
    for line in probands.split('\n'):
        if line:
            key, value = line.split(':')
            probands_dict[key] = value
    return probands_dict

def proband_gt_dict_and_ped_to_inheritance_counts(proband_gt_dict, cohort_gt_dict, ped_file):


    '''
    This script needs to find the line in the ped file for each proband in the dict. 
    First it needs to see if the proband is a trio, and if not it gets added to the unknown inheritance list. 
    If it is a trio, it needs to check if the variant is present in the parents.
    Absence of a variant in the parents is not sufficient, the site must be called with the REF allele in parents for de novo.
    '''

    # Read the PED file which is tab-delimited into a pandas DataFrame

    counts_dict = {"probands": [], "n_monoallelic": 0, "monoallelic_probands": [], "n_biallelic": 0, "biallelic_probands": [], "n_denovo": 0, "denovo_probands" : [], "n_unknown": 0, "unknown_probands": []}

    for key, value in proband_gt_dict.items():
        counts_dict["probands"].append(key)

    ped_df = pd.read_csv(ped_file, sep='\t', header=None, names=['family_id', 'proband_id', 'father_id', 'mother_id', 'sex', 'phenotype'])
    for proband_id, gt in proband_gt_dict.items():
        proband_row = ped_df[ped_df['proband_id'] == proband_id]
        if proband_row.empty:
            #exit with error if proband not found in PED file, since all probands should be in the PED file
            sys.stderr.write(f"Proband {proband_id} not found in PED file.\n")
            counts_dict["n_unknown"] += 1
            counts_dict["unknown_probands"].append(proband_id)
            continue
        else:
            proband_row = proband_row.iloc[0]
        family_id = proband_row['family_id']
        father_id = proband_row['father_id']
        if father_id != "0":
            father_gt = cohort_gt_dict[father_id]
        father_has_variant = father_id in proband_gt_dict.keys()
        mother_id = proband_row['mother_id']
        if mother_id != "0":
            mother_gt = cohort_gt_dict[mother_id]
        mother_has_variant = mother_id in proband_gt_dict.keys()
        # Add to unknown dict if either parent is missing
        if father_id == "0" or mother_id == "0":
            counts_dict["n_unknown"] += 1
            counts_dict["unknown_probands"].append(proband_id)
        elif (father_gt in ["./.", "0", ".", ".|."] or mother_gt in ["./.", "0", ".", ".|."]):
            counts_dict["n_unknown"] += 1
            counts_dict["unknown_probands"].append(proband_id)

        # Only add to the counts if the proband is part of a trio
        elif (father_has_variant or mother_has_variant) and not (father_has_variant and mother_has_variant):
            counts_dict["n_monoallelic"] += 1
            counts_dict["monoallelic_probands"].append(proband_id)
        elif father_has_variant and mother_has_variant:
            counts_dict["n_biallelic"] += 1
            counts_dict["biallelic_probands"].append(proband_id)
        # Check for de novo variants
        elif mother_gt == "0/0" and father_gt == "0/0":
            counts_dict["n_denovo"] += 1
            counts_dict["denovo_probands"].append(proband_id)


    return(counts_dict)

def make_tab_delimited_string(gnomad_id, rest, input_dict):
    """Convert a dictionary into a tab-delimited string with keys in a specified order."""
    
    # Define the desired order of keys
    key_order = [
        "probands", "n_monoallelic", "monoallelic_probands", "n_biallelic", "biallelic_probands", "n_denovo", "denovo_probands", "n_unknown", "unknown_probands"
    ]
    

    # Get the values in the desired order
    ordered_values = [input_dict.get(key, '') for key in key_order]
    
    # Join the values into a tab-delimited string
    return f"{gnomad_id}\t{rest}\t" + "\t".join(map(str, ordered_values))



def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process gnomAD data and convert a semicolon-delimited string into a dictionary.")
    parser.add_argument("gnomad_id", type=str, help="The gnomAD identifier.")
    parser.add_argument("rest", type=str, help="Additional string information.")
    parser.add_argument("probands", type=str, help="A enter delimited string of probands and their GT fields.")
    parser.add_argument("--ped", type=str, help="The path to the PED file.")

    
    # Parse the arguments
    args = parser.parse_args()
    cohort_gt_dict = probands_to_dict(args.probands)
    proband_gt_dict = {key: value for key, value in cohort_gt_dict.items() if value not in ["0/0", "./.", "0|0", "0", ".", ".|."]}
    counts_dict = proband_gt_dict_and_ped_to_inheritance_counts(proband_gt_dict, cohort_gt_dict, args.ped)
    print(make_tab_delimited_string(args.gnomad_id, args.rest, counts_dict))

if __name__ == "__main__":
    main()