#!/usr/bin/env python

import random
from Bio import Phylo
import argparse

def sample_newick(input_file_path, output_file_path, num_samples):
    tree = Phylo.read(input_file_path, 'newick')
    tips = tree.get_terminals()
    sampled_tips = random.sample(tips, int(num_samples))
    sampled_tip_names = [tip.name for tip in sampled_tips]
    for tip in tree.get_terminals():
        if tip.name not in sampled_tip_names:
            tree.prune(tip)
    Phylo.write(tree, output_file_path, 'newick')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tree Sampling. Create a new tree by sampling a subset (count) of tips from the input tree.")
    parser.add_argument('--inp', type=str, help='input file path')
    parser.add_argument('--out', type=str, help='output file path')
    parser.add_argument('--count', type=str, help='number of tips to sample')
    args = parser.parse_args()

    input_file_path = args.inp
    output_file_path = args.out
    num_samples = args.count
    sample_newick(input_file_path, output_file_path, num_samples)