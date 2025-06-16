#!/usr/bin/env python

import random
from Bio import SeqIO
import argparse

def print_random_sequences(fasta_file, num_sequences):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    random_sequences = random.sample(sequences, num_sequences)
    
    for seq in random_sequences:
        print(f">{seq.id}\n{seq.seq}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print num random sequences from a FASTA file.")
    parser.add_argument('--inp', type=str, required=True, help='input FASTA file path')
    parser.add_argument('--num', type=int, required=True, help='number of random sequences to print')
    args = parser.parse_args()

    print_random_sequences(args.inp, args.num)