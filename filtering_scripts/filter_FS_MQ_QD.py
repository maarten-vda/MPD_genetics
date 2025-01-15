import argparse
import pandas as pd

# Function to filter the DataFrame based on given conditions
def filter_dataframe(df):
    # Convert columns to integers, replace errors with NaN, then drop rows where any of these columns are NaN
    for col in ['MQ', 'QD', 'FS', 'phylop']:
        df[col] = pd.to_numeric(df[col], errors='coerce')  # Convert to numeric, coercing invalid values to NaN
    
    # Drop rows where any of the columns 'MQ', 'QD', or 'FS' contain NaN
    df.dropna(subset=['MQ', 'QD', 'FS', 'phylop'], inplace=True)
    
    # Apply the filtering conditions
    filtered_df = df[(df['MQ'] >= 30) & (df['QD'] >= 2) & (df['FS'] <= 60) & (df['phylop'] >= 0) & (df['consequence'] != 'synonymous_variant')]
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
