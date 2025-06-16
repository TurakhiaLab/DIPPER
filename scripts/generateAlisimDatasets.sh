# generateDatasets.sh 2&> generate_100.log
# export PATH=$PATH:/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/scripts
SCRIPTS=/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/scripts
DATASET_SIZE=2000000
LEN=10000
bl=$1
THREADS=32
# DATASET_DIR=/data/swalia/dipper/alisim/datasets/size_50k_length_v
# DATASET_DIR=/data/swalia/dipper/alisim/datasets/bl_effect
DATASET_DIR=/data/swalia/dipper/alisim/datasets/large
# FINAL_DIR=${DATASET_DIR}/${LEN}
# FINAL_DIR=${DATASET_DIR}/len_${LEN}_size_${DATASET_SIZE}_bl_${bl}
FINAL_DIR=${DATASET_DIR}/2M

FASTA_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.unaligned.fa
MSA_FASTA_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.aligned.fa
DIST_PHY_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.dist
DIST_PHY_FILE_MASH=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.mash.dist
DIST_XML_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.xml # fastdist ../dataset/1000/len_1000.fasta > len_1000.xml
MSA_STH_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.stockholm.fa
MSA_PHY_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.phylip.fa
MSA_BACKBONE_FASTA_FILE=$FINAL_DIR/len_${LEN}_leaves_${DATASET_SIZE}.backbone1000.fa
MASH_DIR=$FINAL_DIR/mash


# mkdir -p ${DATASET_DIR}/${DATASET_SIZE}
mkdir -p ${FINAL_DIR}

# 1. Scale Branch lengths using b_scale.py (avg bl=0.1)

# 2. Generate sequences using alisim (source alisim.sh $leaves)
cd ${FINAL_DIR}
/usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim ${FINAL_DIR}/data --length $LEN -m "GTR+G+I" --indel-size POW{1.7/50},POW{1.7/50}  --indel 0.03,0.09 seed 1 -t "RANDOM{yh/${DATASET_SIZE}}" -nt 1 -rlen 0.000001 0.00002 0.00001;
# /usr/bin/time -v timeout 1d /                                                                                                                                                                                                                                    0.000002 0.00002 0.0002;
# timeout 1d /usr/bin/time -f "Alisim Generation Time:\t%e s" ${SCRIPTS}/vmPeak.sh ${SCRIPTS}/alisim.sh $DATASET_SIZE ${FINAL_DIR} $LEN $bl
# mv ${FINAL_DIR}/len_${LEN}_leaves_${DATASET_SIZE}.phy $MSA_PHY_FILE

# # 3. Convert to fasta and stockholm format using format_converter.py
# ${SCRIPTS}/format_converter.py --inp $MSA_PHY_FILE --in_format phylip --out_format fasta --out $MSA_FASTA_FILE
# format_converter.py --inp $MSA_PHY_FILE --in_format phylip --out_format stockholm --out $MSA_STH_FILE

# # 4. split sequences into individual files for mashtree
# mkdir -p $MASH_DIR
# mash_split.py --in_fasta $FASTA_FILE --out-dir $MASH_DIR

# # # 5. Generate distance matrix using mash
# export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
# MASH_TEMP_DIR=mash_temp_dir
# timeout 1d /usr/bin/time -f "Mash distance Matrix generation time:\t%e s" vmPeak.sh /home/swalia@AD.UCSD.EDU/tools/mashtree-1.4.6/bin/mashtree --sketch-size 1000 --kmerlength 15 --numcpus $THREADS --outmatrix $DIST_PHY_FILE_MASH --tempdir $MASH_TEMP_DIR $MASH_DIR/*
# rm -r $MASH_TEMP_DIR

# # 6. Generate distance matrix using DIPPER-Mash
# echo "#Mash-GPU distance matrix"
# timeout 1d /usr/bin/time -f "Mash-GPU distance Matrix generation time:\t%e s" vmPeak.sh /home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/phylo-accel/build/test-mash-gpu -f $FASTA_FILE -s 1000 > $DIST_PHY_FILE

# # 7. Generate distance matrix using mash-triangle
# timeout 1d /usr/bin/time -f "Mash-triangle distance Matrix generation time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh mash triangle -p 32 -k 15 -s 1000  $MASH_DIR/* > mash-tri.dist
# # 7. Generate DIST XML file for fnj
# # conda activate fastphylo
# # timeout 1d /usr/bin/time -f "Mash-GPU distance Matrix generation time:\t%e s" vmPeak.sh fastdist -e -I phylip $DIST_PHY_FILE > $DIST_XML_FILE

# # 8. Split first 1000 sequences for backbone trees for bscampp
# #echo "Splitting first 1000 sequences for backbone trees"
# #bscampp_input.py --inp $MSA_FASTA_FILE --num 1000 > $MSA_BACKBONE_FASTA_FILE
