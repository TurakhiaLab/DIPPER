
LEN=10000
DATASET_SIZE=1000000
THREADS=32
SCRIPTS_DIR=/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/scripts/
# RESULTS_DIR=/data/swalia/dipper/alisim/results/len_${LEN}_leaves_${DATASET_SIZE}_threads_${THREADS}
RESULTS_DIR=/data/swalia/dipper/RNAsim/results/${DATASET_SIZE}_threads_${THREADS}

RAPIDNJ_MSA=$RESULTS_DIR/rapidnj_msa.nwk
FASTME_MSA=$RESULTS_DIR/fastme_msa.nwk
DECENTTREE_MSA=$RESULTS_DIR/decenttree_msa.nwk
QUICKTREE_MSA=$RESULTS_DIR/quicktree_msa.nwk
CCPHYLO_MSA=$RESULTS_DIR/ccphylo_msa.nwk
DIPPER_MSA=$RESULTS_DIR/dipper_msa.nwk
DIPPER_RAW=$RESULTS_DIR/dipper_raw.nwk
DIPPER_RAW_NEW=$RESULTS_DIR/dipper_raw_redo.nwk

RAPIDNJ_DIST=$RESULTS_DIR/rapidnj_dist.nwk
FASTME_DIST=$RESULTS_DIR/fastme_dist.nwk
DECENTTREE_DIST=$RESULTS_DIR/decenttree_dist.nwk
QUICKTREE_DIST=$RESULTS_DIR/quicktree_dist.nwk
CCPHYLO_DIST=$RESULTS_DIR/ccphylo_dist.nwk
DIPPER_DIST=$RESULTS_DIR/dipper_dist.nwk
APPLES_DIST=$RESULTS_DIR/apples.tog.tre

MAPLE=/home/swalia@AD.UCSD.EDU/tools/MAPLE/MAPLEv0.3.6.py
MAPLE_OUT=$RESULTS_DIR/maple_out
CAT_FILE=${MAPLE_OUT}_RFdistances.txt
# REF_TREE=/data/swalia/dipper/alisim/trees/nleaves${DATASET_SIZE}.treefile
REF_TREE=/data/swalia/dipper/RNAsim/dataset/${DATASET_SIZE}/true_tree.nwk
nRF_FILE=$RESULTS_DIR/nRF.txt
rm $nRF_FILE && touch $nRF_FILE

# DecentTree
rm $CAT_FILE
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $DECENTTREE_MSA --output $MAPLE_OUT --overwrite
echo "Decenttree MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $DECENTTREE_DIST --output $MAPLE_OUT --overwrite
# echo "Decenttree DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# FastME
rm $CAT_FILE
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $FASTME_MSA --output $MAPLE_OUT --overwrite
echo "FastME MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $FASTME_DIST --output $MAPLE_OUT --overwrite
# echo "FastME DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# CCPHYLO
rm $CAT_FILE
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $CCPHYLO_MSA --output $MAPLE_OUT --overwrite
echo "CCPHYLO MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $CCPHYLO_DIST --output $MAPLE_OUT --overwrite
# echo "CCPHYLO DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE


# RapidNJ
rm $CAT_FILE
# python3 $SCRIPTS_DIR/parse_rapidnj.py --inp $RAPIDNJ_MSA --out $RESULTS_DIR/rapidnj_msa_trim.nwk --type rapidnj
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/rapidnj_msa_trim.nwk --output $MAPLE_OUT --overwrite
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/rapidnj_msa.nwk --output $MAPLE_OUT --overwrite
echo "RapidNJ MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# rm $CAT_FILE
# python3 $SCRIPTS_DIR/parse_rapidnj.py --inp $RAPIDNJ_DIST --out $RESULTS_DIR/rapidnj_dist_trim.nwk --type rapidnj
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/rapidnj_dist_trim.nwk --output $MAPLE_OUT --overwrite
# echo "RapidNJ DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# QuickTree
rm $CAT_FILE
# python3 $SCRIPTS_DIR/parse_rapidnj.py --inp $QUICKTREE_MSA --out $RESULTS_DIR/quicktree_msa_trim.nwk --type quicktree
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/quicktree_msa_trim.nwk --output $MAPLE_OUT --overwrite
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/quicktree_msa.nwk --output $MAPLE_OUT --overwrite
echo "QuickTree MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# rm $CAT_FILE
# python3 $SCRIPTS_DIR/parse_rapidnj.py --inp $QUICKTREE_DIST --out $RESULTS_DIR/quicktree_dist_trim.nwk --type quicktree
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/quicktree_dist_trim.nwk --output $MAPLE_OUT --overwrite
# echo "QuickTree DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# APPLES
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $APPLES_MSA --output $MAPLE_OUT --overwrite
# echo "APPLES MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $APPLES_DIST --output $MAPLE_OUT --overwrite
# echo "APPLES DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

rm $CAT_FILE
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $RESULTS_DIR/fnj_msa.nwk --output $MAPLE_OUT --overwrite
echo "Fastphylo DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# Dipper
rm $CAT_FILE
python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $DIPPER_MSA --output $MAPLE_OUT --overwrite
echo "Dipper MSA nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# echo "Dipper DIST nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE
# rm $CAT_FILE
# python3 $MAPLE --inputTree $REF_TREE --inputRFtrees $DIPPER_RAW --output $MAPLE_OUT --overwrite
# echo "Dipper RAW nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# rm $CAT_FILE
# python3 $MAPLE --inputTree $DIPPER_RAW_NEW --inputRFtrees $DIPPER_RAW --output $MAPLE_OUT --overwrite
# echo "Dipper RAW NEW nrf: $(awk 'NR==2 {print $2}' $CAT_FILE)" >> $nRF_FILE

# Print nRF values
cat $nRF_FILE
