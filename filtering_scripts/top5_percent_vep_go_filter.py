import argparse
import pandas as pd
import ast

def filter_condition(row):
    #This function just filters all monoallelic records if they are not in X-linked genes
    inheritance_list = ast.literal_eval(row['inheritance'])

    # Check if 'X-linked' or 'X-linked recessive' is NOT in the inheritance list
    if 'X-linked' in inheritance_list or 'X-linked recessive' in inheritance_list:
        return True
    elif row['n_monoallelic'] < 1:
        return True
    else:
        return False  # If the condition is not met, return False

# Function to filter the DataFrame based on given conditions
def filter_dataframe(df):
    # Convert columns to integers, replace errors with NaN, then drop rows where any of these columns are NaN
    for col in ['MQ', 'QD', 'FS', 'phylop']:
        df[col] = pd.to_numeric(df[col], errors='coerce')  # Convert to numeric, coercing invalid values to NaN
    
    # Drop rows where any of the columns 'MQ', 'QD', or 'FS' contain NaN
    df.dropna(subset=['MQ', 'QD', 'FS', 'phylop'], inplace=True)
    
    VEP = "CPT"

    # Apply the filtering conditions
    filtered_df = df[(df[VEP] != ".")]

    # Convert VEP, F_enrichment, C_enrichment, and P_enrichment to numeric, coercing invalid values to NaN
    for col in [VEP, 'F_enrichment', 'C_enrichment', 'P_enrichment']:
        filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')

    top_5_percent_threshold = filtered_df[VEP].quantile(0.95)  ## Swap this depending on VEP more positive or more negative predictions correspond to pathogenic

    filtered_df = filtered_df[(filtered_df[VEP] > top_5_percent_threshold) & (filtered_df['F_enrichment'] > 0.75) & (filtered_df['C_enrichment'] > 0.75) & (filtered_df['P_enrichment'] > 0.75)]
    filtered_df = filtered_df.loc[filtered_df.apply(filter_condition, axis=1)]

    # Sort by VEP in descending order
    filtered_df.sort_values(by=VEP, ascending=False, inplace=True)

    return filtered_df

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Filter rows in a TSV based on specific conditions.')
    parser.add_argument('--input', required=True, help='Input TSV file path')
    parser.add_argument('--output', required=True, help='Output TSV file path')

    args = parser.parse_args()

    # Load the TSV file into a pandas DataFrame
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error loading input file: {e}")
        return

    # Filter the DataFrame based on the conditions
    filtered_df = filter_dataframe(df)

    # Write the filtered DataFrame to the output TSV file
    try:
        filtered_df.to_csv(args.output, sep='\t', index=False)
        print(f"Filtered data saved to {args.output}")
    except Exception as e:
        print(f"Error saving output file: {e}")

if __name__ == '__main__':
    main()
