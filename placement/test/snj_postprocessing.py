from Bio import SeqIO
import argparse
import re



parser = argparse.ArgumentParser(description="dataset_size")
parser.add_argument('--dataset_size', type=int, help='dataset_size')
args = parser.parse_args()

fasta_file = f'/data/zec022/phastsim_datasets/dataset_{args.dataset_size}/sars-cov-2_simulation_output.fasta'

sequence_names = []

for record in SeqIO.parse(fasta_file, "fasta"):
    sequence_names.append(record.id)

# print(sequence_names)

input_file = '/home/zec022@AD.UCSD.EDU/SNJ/results/msa_SNJ_tree_0.nw'
output_file = f'/data/zec022/SNJ/dataset_{args.dataset_size}/tree.nwk'

with open(input_file, 'r') as file:
    content = file.read()

def replace_match(match):
    prefix = match.group(1)  # "(" 或 ","
    index = int(match.group(2))
    if 0 <= index < len(sequence_names):
        return f"{prefix}{sequence_names[index]}"
    return match.group(0)

updated_content = re.sub(r'(\(|,)(\d+)', replace_match, content)

with open(output_file, 'w') as file:
    file.write(updated_content)

# print(updated_content)