import argparse
import csv

def format_proband_ids(input_vcf):
    with open(input_vcf, "r") as f:
        for line in f:
            # Skip lines that are comments
            if line.startswith("##"):
                continue
            # Get the header line and initialize the formatted VCF string
            elif line.startswith("#") and not line.startswith("##"):
                header = line.strip().split("\t")
                formatted_vcf = "\t".join(header[:8]) + "\n"
            else:
                # Process all the lines that are not comments or headers
                line = line.strip().split("\t")
                if line[6] != "PASS":
                    continue
                chrom = line[0]
                chrom_numerical = chrom[3:] if chrom.startswith("chr") else chrom
                pos = line[1]
                ref = line[3]
                alt = line[4]
                genotype_format = line[8].split(":")
                genotype_field = genotype_format.index("GT")
                depth_field = genotype_format.index("DP")
                proband_info = line[9:]
                alleles_proband_dict = {}
                # Make a dictionary containing which proband IDs have which disease allele
                for proband_index, proband in enumerate(proband_info):
                    proband = proband.split(":")
                    # Ignore probands that are homozygous for the reference allele or have missing genotypes
                    if proband[genotype_field] == "./." or proband[0] == "0/0":
                        continue
                    else:
                        genotype = proband[genotype_field]
                        alleles = genotype.split("/")
                        # Ensure alleles has at least two elements
                        if len(alleles) >= 2:
                            allele1 = alleles[0] if alleles[0] != '0' else None
                            allele2 = alleles[1] if alleles[1] != '0' else None
                        else:
                            # Handle unexpected cases where alleles doesn't have at least two elements
                            allele1 = None
                            allele2 = None
                            # Optionally, print a message or log this unexpected case for debugging
                            print(f"Unexpected genotype format: {genotype}")
                        read_depth = proband[depth_field]
                        proband_id = header[proband_index + 9]  # Correct index for proband ID
                        info_field = [proband_id, genotype, read_depth, fetch_diagnosis(proband_id, args.probands)]
                        # Append proband ID to dictionary mapping alleles to proband IDs
                        if allele1 and allele1 != '0':
                            if allele1 not in alleles_proband_dict:
                                alleles_proband_dict[allele1] = []
                            alleles_proband_dict[allele1].append(info_field)
                        if allele2 and allele2 != '0':
                            if allele2 not in alleles_proband_dict:
                                alleles_proband_dict[allele2] = []
                            alleles_proband_dict[allele2].append(info_field)
                # Check if the alt allele has multiple values, if it does split them
                if "," in alt:
                    alts = alt.split(",")
                else:
                    alts = [alt]
                for allele_idx, alt_allele in enumerate(alts):
                    # Generate the info field for the current alt allele
                    minimal_info = line[:9]
                    if alt_allele != "*":
                        minimal_info[2] = f"{chrom_numerical}-{pos}-{ref}-{alt_allele}"  # Set the ID to gnomAD style ID
                    elif alt_allele == "*":
                        minimal_info[2] = f"{chrom_numerical}-{pos}-{ref}-del"
                    minimal_info[4] = alt_allele
                    minimal_info[7] = "samples=" + str(alleles_proband_dict[str(allele_idx+1)])
                    ##remove the format field
                    minimal_info.pop()
                    line_string = "\t".join(minimal_info)
                    formatted_vcf += line_string + "\n"
    return formatted_vcf
                    
##Function to fetch whether a proband is diagnosed or undiagnosed
def fetch_diagnosis(proband_id, proband_file):
    # Convert search_string to lowercase for case-insensitive search
    search_string = proband_id.lower()
    
    # Open the CSV file
    with open(proband_file, 'r', newline='') as file:
        reader = csv.reader(file)

        # Iterate through each row in the CSV
        for row in reader:
            # Convert each cell to lowercase for case-insensitive comparison
            row_lower = [cell.lower() for cell in row]

            # Iterate through each cell in the row
            for col_index, cell in enumerate(row_lower):
                # Check if the search_string exists in the cell
                if search_string in cell:
                    # Determine the value based on the column index
                    if col_index == 0:
                        # If the search string is in the first column
                        return row[3]
                    elif col_index == 1:
                        # If the search string is in the second column
                        return f"father_{row[3]}"
                    elif col_index == 2:
                        # If the search string is in the third column
                        return f"mother_{row[3]}"
    # If the search string is not found
    return "."

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve which proband IDs have the disease allele and format them")
    parser.add_argument("--input", required=True, help="Path to input VCF file")
    parser.add_argument("--probands", required=True, help="Path to proband IDs file")
    args = parser.parse_args()
    print(format_proband_ids(args.input))