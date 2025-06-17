HOME_DIR="/home/swalia@AD.UCSD.EDU/work/phylo-accel/build"
DATA_DIR="/home/swalia@AD.UCSD.EDU/work/phylo-accel/build"
cd $HOME_DIR

# Make directories
mkdir -p results/distMat
mkdir -p results/newick

fasta=$1

time ./test-mash-gpu -e 20 -f $fasta -s 10000 > results/distMat/drosophila_20.phy
time ./test-nj -i results/distMat/drosophila_20.phy > results/newick/drosophila_20.nwk
    
