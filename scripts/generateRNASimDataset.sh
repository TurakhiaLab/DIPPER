# generateDatasets.sh 2&> generate_100.log
export PATH=$PATH:~/work/phylo-accel/scripts
DATASET_SIZE=200000
THREADS=32
DATASET_DIR=/data/swalia/dipper/RNAsim/dataset

MSA_FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.fa
MSA_STH_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.stockholm.fa
MSA_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.phylip.fa

format_converter.py --inp $MSA_FASTA_FILE --in_format fasta --out_format phylip --out $MSA_PHY_FILE
format_converter.py --inp $MSA_FASTA_FILE --in_format fasta --out_format stockholm --out $MSA_STH_FILE

# mkdir -p $MASH_DIR
# mash_split.py --in_fasta $FASTA_FILE --out-dir $MASH_DIR

# # Generate distance matrix using mash
# export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
# MASH_TEMP_DIR=mash_temp_dir
# timeout 1d /usr/bin/time -f "Mash distance Matrix generation time:\t%e s" vmPeak.sh /home/swalia@AD.UCSD.EDU/tools/mashtree-1.4.6/bin/mashtree --sketch-size 1000 --kmerlength 15 --numcpus $THREADS --outmatrix $DIST_PHY_FILE_MASH --tempdir $MASH_TEMP_DIR $MASH_DIR/*
# rm -r $MASH_TEMP_DIR

# echo "#Mash-GPU distance matrix"
# timeout 1d /usr/bin/time -f "Mash-GPU distance Matrix generation time:\t%e s" vmPeak.sh /home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/phylo-accel/build/test-mash-gpu -f $FASTA_FILE -s 1000 > $DIST_PHY_FILE

# #  Split first 1000 sequences for backbone trees for bscampp
# echo "Splitting first 1000 sequences for backbone trees"
# bscampp_input.py --inp $MSA_FASTA_FILE --num 1000 > $MSA_BACKBONE_FASTA_FILE
