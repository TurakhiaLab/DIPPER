from Bio import SeqIO
from sys import argv

import pandas as pd

def remove_queries_from_matrix(matrix_file, remove_ids_file, output_file):
    # Read the list of query IDs to remove
    print(f"Reading IDs to remove from {remove_ids_file}")
    with open(remove_ids_file, 'r') as f:
        remove_ids = set(line.strip() for line in f if line.strip())

    print(f"Found {len(remove_ids)} IDs to remove.")
    # Load the matrix as a DataFrame
    df = pd.read_csv(matrix_file, sep='\t', index_col=0)

    print(f"Loaded matrix with shape {df.shape}.")
    # Remove the rows that match the query IDs
    filtered_df = df[~df.index.isin(remove_ids)]

    print(f"Filtered matrix shape: {filtered_df.shape}.")
    # Write the filtered matrix to a new file
    filtered_df.to_csv(output_file, sep='\t')

# Example usage
matrix_file = argv[1]  
remove_ids_file = argv[2]
output_file = argv[3]
if (len(argv) != 4):
    print("Usage: python remove_queries_from_matrix.py <matrix_file> <remove_ids_file> <output_file>")

remove_queries_from_matrix(matrix_file, remove_ids_file, output_file)