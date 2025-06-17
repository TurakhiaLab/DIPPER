from Bio import SeqIO
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument('--input_msa', type=str, help='file path to the input msa')
parser.add_argument('--output', type=str, help='file path to the output')
args = parser.parse_args()

fasta_file = args.input_msa
sequences = []

for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(list(str(record.seq)))

data = np.array(sequences)

np.save(args.output, data)