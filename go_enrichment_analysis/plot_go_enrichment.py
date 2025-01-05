import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
# Use Agg backend
matplotlib.use('Agg')

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Generate bar chart from GO terms with enrichment > 20.")
    parser.add_argument("--p", required=True, help="Path to TSV file for PANTHER GO-Slim Biological Process")
    parser.add_argument("--c", required=True, help="Path to TSV file for PANTHER GO-Slim Cellular Component")
    parser.add_argument("--f", required=True, help="Path to TSV file for PANTHER GO-Slim Molecular Function")
    parser.add_argument("--output", required=True, help="Path to save the output bar chart")
    args = parser.parse_args()

    # Load dataframes
    p_df = pd.read_csv(args.p, sep='\t')
    c_df = pd.read_csv(args.c, sep='\t')
    f_df = pd.read_csv(args.f, sep='\t')

    # Convert fold enrichment to numeric
    p_df['upload_1 (fold Enrichment)'] = pd.to_numeric(p_df['upload_1 (fold Enrichment)'], errors='coerce')
    c_df['upload_1 (fold Enrichment)'] = pd.to_numeric(c_df['upload_1 (fold Enrichment)'], errors='coerce')
    f_df['upload_1 (fold Enrichment)'] = pd.to_numeric(f_df['upload_1 (fold Enrichment)'], errors='coerce')

    # Filter dataframes 
    p_filtered = p_df[
        (p_df['upload_1 (fold Enrichment)'] > 20) & 
        (p_df['Homo sapiens - REFLIST (20580)'] >= 20) & 
        (p_df['Homo sapiens - REFLIST (20580)'] <= 200) & 
        (p_df['upload_1 (414)'] >= 10)
    ]

    # Filter dataframes for C
    c_filtered = c_df[
        (c_df['Homo sapiens - REFLIST (20580)'] >= 7) & 
        (c_df['Homo sapiens - REFLIST (20580)'] <= 70) & 
        (c_df['upload_1 (fold Enrichment)'] > 20)
    ]

    # Filter dataframes for F
    f_filtered = f_df[
        (f_df['upload_1 (fold Enrichment)'] > 20) & 
        (f_df['Homo sapiens - REFLIST (20580)'] >= 5) & 
        (f_df['Homo sapiens - REFLIST (20580)'] <= 40)
    ]

    # This suppresses the error when using a slice of a dataframe
    p_filtered = p_filtered.copy()
    c_filtered = c_filtered.copy()
    f_filtered = f_filtered.copy()

    # Sort by the fold enrichment for each of the subsets
    p_filtered.sort_values(by='upload_1 (fold Enrichment)', ascending=True, inplace=True)
    c_filtered.sort_values(by='upload_1 (fold Enrichment)', ascending=True, inplace=True)
    f_filtered.sort_values(by='upload_1 (fold Enrichment)', ascending=True, inplace=True)

    # Plot setup
    fig, ax = plt.subplots(figsize=(10, 8))


    # Adjust y-positions for stacking categories without overlap
    p_positions = range(len(p_filtered))
    c_positions = range(len(p_filtered), len(p_filtered) + len(c_filtered))
    f_positions = range(len(p_filtered) + len(c_filtered), len(p_filtered) + len(c_filtered) + len(f_filtered))

    # Generate horizontal bar plot
    ax.barh(p_positions, p_filtered['upload_1 (fold Enrichment)'], color='blue', label='P', height=0.4)
    ax.barh(c_positions, c_filtered['upload_1 (fold Enrichment)'], color='green', label='C', height=0.4)
    ax.barh(f_positions, f_filtered['upload_1 (fold Enrichment)'], color='red', label='F', height=0.4)

    # Add horizontal lines between categories
    ax.axhline(y=max(p_positions) + 0.5, color='black', xmin=-10, xmax=1, linewidth=0.8, linestyle='-')
    ax.axhline(y=max(c_positions) + 0.5, color='black', xmin=-10, xmax=1, linewidth=0.8, linestyle='-')


    # Set y-ticks and labels
    all_labels = list(p_filtered.iloc[:, 0]) + list(c_filtered.iloc[:, 0]) + list(f_filtered.iloc[:, 0])
    ax.set_yticks(list(p_positions) + list(c_positions) + list(f_positions))
    ax.set_yticklabels(all_labels)

    # Add labels and legend
    ax.set_xlabel('Fold Enrichment')
    ax.set_ylabel('GO Term')
    ax.legend(title="Category", loc="lower right")

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(args.output)
    print(f"Bar chart saved to {args.output}")

if __name__ == "__main__":
    main()
