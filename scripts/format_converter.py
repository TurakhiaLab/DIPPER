#!/usr/bin/env python
from Bio import SeqIO
from Bio import AlignIO
import argparse

def convert_file(input_file_path, output_file_path, input_format, output_format, chunk_size=1000):
    processed=0
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        records = SeqIO.parse(input_file, input_format)
        batch = []
        for record in records:
            batch.append(record)
            processed += 1
            if len(batch) >= chunk_size:
                SeqIO.write(batch, output_file, output_format)
                batch = []
            # print(processed, batch.__len__(), chunk_size)
        if batch:
            SeqIO.write(batch, output_file, output_format)
        
def phylip_to_fasta(input_file_path, output_file_path, batch_size=20000):
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        first_line = input_file.readline().strip()
        num_sequences, sequence_length = map(int, first_line.split())

        sequences = {}
        processed = 0   
        for _ in range(num_sequences):
            line = input_file.readline().strip()
            parts = line.split()
            identifier = parts[0]
            sequence = ''.join(parts[1:])
            sequences[identifier] = sequence
            processed += 1
            if (processed == batch_size):
                for identifier, sequence in sequences.items():
                    output_file.write(f'>{identifier}\n')
                    output_file.write(sequence + '\n')
                sequences = {}
                processed = 0
        if sequences:
            for identifier, sequence in sequences.items():
                output_file.write(f'>{identifier}\n')
                output_file.write(sequence + '\n')


parser = argparse.ArgumentParser(description="")
parser.add_argument('--inp', type=str, help='input file path')
parser.add_argument('--in_format', type=str, help='input file format (stockholm/phylip/fasta)')
parser.add_argument('--out', type=str, help='output file path')
parser.add_argument('--out_format', type=str, help='output file format (stockholm/phylip/fasta)')
args = parser.parse_args()

if ((args.in_format != "stockholm") & (args.in_format != "fasta") & (args.in_format != "phylip")):
    print("Input format error")
    print(args.in_format)
    print(args.out_format)
    exit(0)
if ((args.out_format != "stockholm") & (args.out_format != "fasta") & (args.out_format != "phylip")):
    print("Output format error")
    exit(0)
if (args.in_format == args.out_format):
    exit(0)
if ((args.in_format == "") or (args.out_format == "")):
    print("Input/output format error")
    exit(0)
records = SeqIO.parse(args.inp, args.in_format)
SeqIO.write(records, args.out, args.out_format)
# # phylip_to_fasta(args.inp, args.out)
# convert_file(args.inp, args.out, args.in_format, args.out_format)