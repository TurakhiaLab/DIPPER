from Bio import Phylo

def bl_scale(input_file, output_file, scale_factor):
    tree = Phylo.read(input_file, 'newick')
    
    for clade in tree.find_clades():
        if clade.branch_length:
            clade.branch_length /= scale_factor
            if (clade.branch_length <= 0):
                clade.branch_length = 0.0000001

            #elif (clade.branch_length < 0.0001):
             #   clade.branch_length *= scale_factor
    
    Phylo.write(tree, output_file, 'newick')

def bl_range(input_file):
    tree = Phylo.read(input_file, 'newick')
    branch_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length not in (None, 0)]
    print("No. of tips:", len(tree.get_terminals()))
    if branch_lengths:
        min_length = min(branch_lengths)
        max_length = max(branch_lengths)
        avg_length = sum(branch_lengths) / len(branch_lengths)
        print(f"Branch length range: {min_length} - {max_length}")
        print(f"Avg Branch length : {avg_length}")
    else:
        print("No branch lengths found in the tree.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Scale branch lengths in a Newick file.")
    parser.add_argument("--inp", help="Path to the input Newick file")
    parser.add_argument("--out", help="Path to the output Newick file")
    parser.add_argument("--scale", type=float, help="Factor by which to unscale branch lengths")
    parser.add_argument("--range", type=bool, help="Print the range of branch lengths")

    args = parser.parse_args()
    if (args.range):
        bl_range(args.inp)
        exit(0)
    bl_scale(args.inp, args.out, args.scale)
