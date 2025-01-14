from Bio import SeqIO
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="dataset_size")
parser.add_argument('--dataset_size', type=int, help='dataset_size')
args = parser.parse_args()

fasta_file = f'/data/zec022/phastsim_datasets/dataset_{args.dataset_size}/sars-cov-2_simulation_output.fasta'
sequences = []

for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(list(str(record.seq)))

data = np.array(sequences)

np.save(f'/home/zec022@AD.UCSD.EDU/SNJ/data/msa.npy', data)
