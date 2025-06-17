#!/bin/bash

# Command to use Astral Pro 2: /home/swalia@AD.UCSD.EDU/tools/ASTER-Linux/bin/astral-pro2 -t 8 -u 3 -i gene_tree_merged.nwk -o output.nwk -a mapping.txt

HOME_DIR="/home/swalia@AD.UCSD.EDU/work/phylo-accel/build"
DATA_DIR="/home/swalia@AD.UCSD.EDU/work/phylo-accel/build"
cd $HOME_DIR

# Make directories
mkdir -p results/distMat
mkdir -p results/newick

# Set Seeds
SEEDS=""
SEEDS+="20 "
SEEDS+="40 "
SEEDS+="60 "
# SEEDS+="80 "
# SEEDS+="100 "
# SEEDS+="120 "
# SEEDS+="140 "
# SEEDS+="160 "
# SEEDS+="180 "
# SEEDS+="200 "

fasta=$1

# Step 1: Generate mutiple trees by modifying hash function
for SEED in $SEEDS
do 
    ./test-mash-gpu -e $SEED -f $fasta > results/distMat/drosophila_20_$SEED.phy
    ./test-nj -i results/distMat/drosophila_20_$SEED.phy > results/newick/drosophila_20_$SEED.nwk
done

# # Step 2: Combine all trees
rm -f results/newick/drosophila_20_combined.nwk
for SEED in $SEEDS
do 
    cat results/newick/drosophila_20_$SEED.nwk >> results/newick/drosophila_20_combined.nwk
done

# /home/swalia@AD.UCSD.EDU/tools/ASTER-Linux/bin/astral-pro2 -t 8 -u 3 -i results/newick/sars_combined.nwk -o output.nwk



    
