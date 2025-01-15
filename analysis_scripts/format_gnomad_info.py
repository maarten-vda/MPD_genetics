import argparse
import sys

def parse_string_to_dict(input_string):
    """Convert a semicolon-delimited string into a dictionary."""
    result = {}
    for pair in input_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            if key.strip().startswith('AF'):
                result[key.strip()] = value.strip()
            if key.strip() == "inbreeding_coeff":
                result[key.strip()] = value.strip()
            if key.strip() == "spliceai_ds_max":
                result[key.strip()] = value.strip()
            if key.strip() == "FS":
                result[key.strip()] = value.strip()
            if key.strip() == "MQ":
                result[key.strip()] = value.strip()
            if key.strip() == "QD":
                result[key.strip()] = value.strip()
            if key.strip() == "phylop":
                result[key.strip()] = value.strip()
            

    if "AF_grpmax" not in result.keys():
        result["AF_grpmax"] = max(value for key, value in result.items() if key.startswith("AF"))

    if "spliceai_ds_max" not in result.keys():
        result["spliceai_ds_max"] = "."

    return result

def make_tab_delimited_string(gnomad_id, rest, input_dict):
    """Convert a dictionary into a tab-delimited string with keys in a specified order."""
    
    # Define the desired order of keys
    key_order = [
        "AF", "AF_XX", "AF_XY", "AF_afr_XX", "AF_afr_XY", "AF_afr", 
        "AF_ami_XX", "AF_ami_XY", "AF_ami", "AF_amr_XX", "AF_amr_XY", "AF_amr", 
        "AF_asj_XX", "AF_asj_XY", "AF_asj", "AF_eas_XX", "AF_eas_XY", "AF_eas", 
        "AF_fin_XX", "AF_fin_XY", "AF_fin", "AF_mid_XX", "AF_mid_XY", "AF_mid", 
        "AF_nfe_XX", "AF_nfe_XY", "AF_nfe", "AF_raw", "AF_remaining_XX", 
        "AF_remaining_XY", "AF_remaining", "AF_sas_XX", "AF_sas_XY", "AF_sas", 
        "AF_grpmax", "FS", "MQ", "QD", "inbreeding_coeff", "spliceai_ds_max", "phylop"
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
    parser.add_argument("gnomad_query", type=str, help="The semicolon-delimited key=value string to parse.")

    # Parse command line arguments
    args = parser.parse_args()

    if args.gnomad_query == "":
        gnomad_dict = {'AF': '.', 'AF_XX': '.', 'AF_XY': '.', 'AF_afr_XX': '.', 'AF_afr_XY': '.', 'AF_afr': '.', 'AF_ami_XX': '.', 'AF_ami_XY': '.', 'AF_ami': '.', 'AF_amr_XX': '.', 'AF_amr_XY': '.', 'AF_amr': '.', 'AF_asj_XX': '.', 'AF_asj_XY': '.', 'AF_asj': '.', 'AF_eas_XX': '.', 'AF_eas_XY': '.', 'AF_eas': '.', 'AF_fin_XX': '.', 'AF_fin_XY': '.', 'AF_fin': '.', 'AF_mid_XX': '.', 'AF_mid_XY': '.', 'AF_mid': '.', 'AF_nfe_XX': '.', 'AF_nfe_XY': '.', 'AF_nfe': '.', 'AF_raw': '.', 'AF_remaining_XX': '.', 'AF_remaining_XY': '.', 'AF_remaining': '.', 'AF_sas_XX': '.', 'AF_sas_XY': '.', 'AF_sas': '.', 'AF_grpmax': '.', 'FS': '.', 'MQ': '.', 'QD': '.', 'inbreeding_coeff': '.', 'spliceai_ds_max': '.', 'phylop': '.'}
        print(make_tab_delimited_string(args.gnomad_id, args.rest, gnomad_dict))
        sys.exit(0)
    else:
        # Parse the semicolon-delimited string into a dictionary
        gnomad_dict = parse_string_to_dict(args.gnomad_query)
        print(make_tab_delimited_string(args.gnomad_id, args.rest, gnomad_dict))
        sys.exit(0)

if __name__ == "__main__":
    main()