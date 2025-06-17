#!/bin/bash
# Ex: source mash.sh mash len_1000.dist
mash_directory=$1
files=("$mash_directory"/*.fa)
num_files=${#files[@]}
kmer=15
THREADS=32

for ((i=0; i<num_files; i++)); do
    a=${files[i]}
    base_a=$(basename "$a" .fa)
    mash sketch -p $THREADS "$a" -o "$base_a" -k $kmer
done

avg_dist=0
a=${files[0]}
base_a=$(basename "$a" .fa)
for ((j=0; j<num_files; j++)); do
    b=${files[j]}
    base_b=$(basename "$b" .fa)
    echo "mash dist -p $THREADS "$base_a.msh" "$base_b.msh""
    dist=$(mash dist -p $THREADS "$base_a.msh" "$base_b.msh" | cut -f3)    
    # avg_dist=$(echo "$avg_dist + $dist" | bc)
    avg_dist=$(awk "BEGIN {print $avg_dist + $dist}")
    echo $avg_dist
done
avg_dist=$(awk "BEGIN {print $avg_dist / $num_files}")
echo "Average distance: $avg_dist"
# rm *.msh

