#!/bin/bash

# TOOLS_DIR=/home/swalia/.conda/envs/dipper_baseline/bin/
# SCRIPTS_DIR=/home/swalia/dipper/scripts
# DATASET_DIR=/home/swalia/dipper/datasets

TOOLS_DIR=/home/swalia@AD.UCSD.EDU/.conda/envs/dipper_baselines/bin
SCRIPTS_DIR=/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/scripts
# DATASET_DIR=/data/swalia/dipper/alisim/datasets
DATASET_DIR=/data/swalia/dipper/RNAsim/dataset
BSCAMPP_DIR=/home/swalia@AD.UCSD.EDU/tools/BSCAMPP-1.0.0/BSCAMPP_code

DATASET_SIZE=$2
LEN=10000
THREADS=32
TIME_UTIL=/usr/bin/time
VMPEAK=$SCRIPTS_DIR/vmPeak.sh

# if dataset size is greater than 100k, ser DIPPER to dipper-dc
DIPPER=/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/build/dipper
if [ $DATASET_SIZE -gt 1000000 ]; then
    DIPPER=/home/swalia@AD.UCSD.EDU/work/phylo-accel-zexing-dev/DIPPER_sumit/DIPPER/build/dipper-dc
fi
MASH_GPU=$TOOLS_DIR/test-mash-gpu
MASHTREE=$TOOLS_DIR/mashtree
FASTTREE=$TOOLS_DIR/FastTree
QUICKTREE=$TOOLS_DIR/quicktree
DECENTTREE=$TOOLS_DIR/decenttree
FASTME=$TOOLS_DIR/fastme
RAPIDNJ=$TOOLS_DIR/rapidnj
FNJ=$TOOLS_DIR/fnj
FASTDIST=$TOOLS_DIR/fastdist
CCPHYLO=$TOOLS_DIR/ccphylo
SNJ=$TOOLS_DIR/snj


# Alisim
# FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.unaligned.fa
# MSA_FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.aligned.fa
# # DIST_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.dist
# DIST_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/attotree.dist
# DIST_XML_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.xml # fastdist ../dataset/1000/len_1000.fasta > len_1000.xml
# MSA_STH_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.stockholm.fa
# MSA_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.phylip.fa
# MSA_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.phylip.fa
# MSA_BACKBONE_FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/len_${LEN}_leaves_${DATASET_SIZE}.backbone1000.fa
# MASH_DIR=$DATASET_DIR/$DATASET_SIZE/mash
# RESULT_DIR=/data/swalia/dipper/alisim/results/len_${LEN}_leaves_${DATASET_SIZE}_threads_${THREADS}

# RNAsim
MSA_FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.fa
MSA_STH_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.stockholm.fa
MSA_PHY_FILE=$DATASET_DIR/$DATASET_SIZE/true_aligned.phylip.fa
MSA_BACKBONE_FASTA_FILE=$DATASET_DIR/$DATASET_SIZE/ture_aligned.backbone1000.fa
MSA_NPY=$DATASET_DIR/$DATASET_SIZE/data/true_aligned.npy
RESULT_DIR=/data/swalia/dipper/RNAsim/results/${DATASET_SIZE}_threads_${THREADS}

mkdir -p $RESULT_DIR

run_type=$1
############# Generate dataset #############
# # Build dist phy file
# echo "Generating distance file in Phylip format..."
# $TIME_UTIL -f "Elapsed time for generating distance file in Phylip format: %E " $MASH_GPU -f $FASTA_FILE > $DIST_PHY_FILE

# # For mashtree, split the sequences into separate files
# echo "Generating dataset for mashtree..."
# mkdir -p $MASH_DIR
# python3 $SCRIPTS_DIR/mash_split.py --in_fasta $FASTA_FILE --out $MASH_DIR

# # Stockholm/Phylip MSA format for quicktree and rapidnj
# echo "Generating Stockholm/Phylip MSA files..."
# $TIME_UTIL -f  "Elapsed time for generating Phylip MSA files: %E " python3 $SCRIPTS_DIR/format_converter.py --in_msa $MSA_FASTA_FILE --out $MSA_PHY_FILE --out_format p
# $TIME_UTIL -f "Elapsed time for generating Stockholm MSA files: %E " python3 $SCRIPTS_DIR/format_converter.py --in_msa $MSA_FASTA_FILE --out $MSA_STH_FILE --out_format s
##################################################################
# Run the Dipper
# echo "Generating Newick tree file using DIPPER..."
# $TIME_UTIL -f "Elapsed time for generating Newick tree file using quicktree from MSA: %E " $QUICKTREE -in a $MSA_STH_FILE > quicktree_msa.nwk
# $TIME_UTIL -f "Elapsed time for generating Newick tree file using quicktree from DIST: %E " $QUICKTREE -in m $DIST_PHY_FILE > quicktree_dist.nwk


if [ $run_type == "msa" ]; then
    # echo "Generating Newick tree file using rapidnj..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $RAPIDNJ $MSA_STH_FILE -i sth -c $THREADS -a jc > ${RESULT_DIR}/rapidnj_msa.nwk
    # python3 $SCRIPTS_DIR/parse_rapidnj.py --inp ${RESULT_DIR}/rapidnj_msa.nwk --out ${RESULT_DIR}/rapidnj_msa_trim.nwk --type rapidnj
    # mv ${RESULT_DIR}/rapidnj_msa.nwk ${RESULT_DIR}/rapidnj_msa.nwk.orig
    # mv ${RESULT_DIR}/rapidnj_msa_trim.nwk ${RESULT_DIR}/rapidnj_msa.nwk


    # echo "Generating Newick tree file using ccphylo..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $CCPHYLO dist -i $MSA_FASTA_FILE -o ${RESULT_DIR}/ccphylo_dist.phy -t $THREADS -C 1 &&  $CCPHYLO tree -i ${RESULT_DIR}/ccphylo_dist.phy -o ${RESULT_DIR}/ccphylo_msa.nwk -t $THREADS
    # rm ${RESULT_DIR}/ccphylo_dist.phy
    
    # echo "Generating Newick tree file using decenttree..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DECENTTREE -fasta $MSA_FASTA_FILE -nt $THREADS -out ${RESULT_DIR}/decenttree_msa.nwk -t NJ-V -no-banner #Default JC model

    echo "Generating Newick tree file using fastme..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK  $FASTME -i $MSA_PHY_FILE -o ${RESULT_DIR}/fastme_msa.nwk -m N -T $THREADS -d J -w n

    echo "Generating Newick tree file using quicktree..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $QUICKTREE -in a $MSA_STH_FILE > ${RESULT_DIR}/quicktree_msa.nwk
    python3 $SCRIPTS_DIR/parse_rapidnj.py --inp ${RESULT_DIR}/quicktree_msa.nwk --out ${RESULT_DIR}/quicktree_msa_trim.nwk --type quicktree
    mv ${RESULT_DIR}/quicktree_msa.nwk ${RESULT_DIR}/quicktree_msa.nwk.orig
    mv ${RESULT_DIR}/quicktree_msa_trim.nwk ${RESULT_DIR}/quicktree_msa.nwk

    # echo "Generating Newick tree file using bscampp..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DIPPER -f $MSA_BACKBONE_FASTA_FILE -a 1 -i m -o t > ${RESULT_DIR}/dipper_backbone1000.nwk
    # cd $BSCAMPP_DIR
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK EPA-ng-BSCAMPP.py -d out -a $MSA_FASTA_FILE  -m JC -t ${RESULT_DIR}/dipper_backbone1000.nwk -i $BSCAMPP_DIR/testing/jc_model
    
    # echo "Generating Newick tree file using fnj..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $FASTDIST $MSA_FASTA_FILE -D JC | $FNJ -O newick > ${RESULT_DIR}/fnj_msa.nwk 
fi 

if [  $run_type == "dist" ]; then
    
    # Run the CCPhylo
    echo "Generating Newick tree file using ccphylo..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $CCPHYLO tree -i $DIST_PHY_FILE -o ${RESULT_DIR}/ccphylo_dist.nwk -t $THREADS 

    # Run the RapidNJ
    echo "Generating Newick tree file using rapidnj..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $RAPIDNJ $DIST_PHY_FILE -i pd -c $THREADS > ${RESULT_DIR}/rapidnj_dist.nwk

    # Run the DecentTree (Many options available)
    echo "Generating Newick tree file using decenttree..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DECENTTREE -in $DIST_PHY_FILE -nt $THREADS -out ${RESULT_DIR}/decenttree_dist.nwk -t NJ-V -no-banner 

    
    # # Run the FastME
    # echo "Generating Newick tree file using fastme..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $FASTME -i $DIST_PHY_FILE -o ${RESULT_DIR}/fastme_dist.nwk -m N -T $THREADS
    
    # # Run the QuickTree
    # echo "Generating Newick tree file using quicktree..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $QUICKTREE -in m $DIST_PHY_FILE > ${RESULT_DIR}/quicktree_dist.nwk


    # echo "Generating Newick tree file using dipper from Raw sequences..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DIPPER -f $FASTA_FILE -a 0 -i r -o t -s 1000 > ${RESULT_DIR}/dipper_raw.nwk
    # Run the FNJ
    # # XML distance format for fnj
    # echo "Generating XML distance file from Phylip format..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK  $FASTDIST -I phylip $DIST_PHY_FILE > $DIST_XML_FILE
    # $TIME_UTIL -f "Elapsed time for generating Newick tree file using fnj from DIST: %E " $FNJ $DIST_XML_FILE -O newick > fnj_dist.nwk

fi 

if [ $run_type == "mash" ]; then
    # Run MashTree
    echo "Generating Newick tree file using mashtree..."
    export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $MASHTREE ---numcpus $THREADS $MASH_DIR/* > mashtree_msa.nwk
fi

if [ $run_type == "dipper" ]; then
    # Run dipper all three modes
    echo "Generating Newick tree file using dipper from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DIPPER -f $MSA_FASTA_FILE -i m -o t -s 1000 -t 2 > ${RESULT_DIR}/dipper_msa.nwk
    #echo "Generating Newick tree file using dipper from Distance Matrix..."
    #timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DIPPER -f $DIST_PHY_FILE -a 0 -i d -o t -s 1000 > ${RESULT_DIR}/dipper_dist.nwk
    # echo "Generating Newick tree file using dipper from Raw sequences..."
    # timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DIPPER -f $FASTA_FILE -a 0 -i r -o t -s 1000 > ${RESULT_DIR}/dipper_raw.nwk
fi

if [ $run_type == "quicktree" ]; then
    # Run quicktree
    echo "Generating Newick tree file using quicktree from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $QUICKTREE -in a $MSA_STH_FILE > ${RESULT_DIR}/quicktree_msa.nwk
    echo "Generating Newick tree file using quicktree from Distance Matrix..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $QUICKTREE -in m $DIST_PHY_FILE > ${RESULT_DIR}/quicktree_dist.nwk
fi

if [ $run_type == "decenttree" ]; then
    # Run decenttree
    echo "Generating Newick tree file using decenttree from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DECENTTREE -fasta $MSA_FASTA_FILE -nt $THREADS -out ${RESULT_DIR}/decenttree_msa.nwk -t NJ-V -no-banner #Default JC model
    echo "Generating Newick tree file using decenttree from Distance Matrix..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $DECENTTREE -in $DIST_PHY_FILE -nt $THREADS -out ${RESULT_DIR}/decenttree_dist.nwk -t NJ-V -no-banner
fi

if [ $run_type == "fastme" ]; then
    # Run fastme
    echo "Generating Newick tree file using fastme from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK  $FASTME -i $MSA_PHY_FILE -o ${RESULT_DIR}/fastme_msa.nwk -m N -T $THREADS -d J -w n
    echo "Generating Newick tree file using fastme from Distance Matrix..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $FASTME -i $DIST_PHY_FILE -o ${RESULT_DIR}/fastme_dist.nwk -m N -T $THREADS -d J -w n
fi

if [ $run_type == "rapidnj" ]; then
    # Run rapidnj
    echo "Generating Newick tree file using rapidnj from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $RAPIDNJ $MSA_STH_FILE -i sth -c $THREADS > ${RESULT_DIR}/rapidnj_msa.nwk
    echo "Generating Newick tree file using rapidnj from Distance Matrix..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $RAPIDNJ $DIST_PHY_FILE -i pd -c $THREADS > ${RESULT_DIR}/rapidnj_dist.nwk
fi

if [ $run_type == "ccphylo" ]; then
    # Run ccphylo
    echo "Generating Newick tree file using ccphylo from MSA..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $CCPHYLO dist -i $MSA_FASTA_FILE -o ${RESULT_DIR}/ccphylo_dist.phy -t $THREADS -C 1 &&  $CCPHYLO tree -i ${RESULT_DIR}/ccphylo_dist.phy -o ${RESULT_DIR}/ccphylo_msa.nwk -t $THREADS
    rm ${RESULT_DIR}/ccphylo_dist.phy
    echo "Generating Newick tree file using ccphylo from Distance Matrix..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $CCPHYLO dist -i $DIST_PHY_FILE -o ${RESULT_DIR}/ccphylo_dist.phy -t $THREADS -C 1 &&  $CCPHYLO tree -i ${RESULT_DIR}/ccphylo_dist.phy -o ${RESULT_DIR}/ccphylo_dist.nwk -t $THREADS
    rm ${RESULT_DIR}/ccphylo_dist.phy
fi

if [ $run_type == "snj" ]; then
    # Run ccphylo
    echo "Generating Newick tree file using SNJ from MSA..."
    # mkdir /data/zec022/SNJ/dataset_$n
    mkdir -p $DATASET_DIR/$DATASET_SIZE/data
    python3.8 ${SCRIPTS_DIR}/snj_msa_input.py --input_msa $MSA_FASTA_FILE --output ${MSA_NPY}
    log_n=$(echo "l($DATASET_SIZE)" | bc -l)
    floor_log_n=$(echo "$log_n / 1" | bc)
    result=$(echo "scale=10; sqrt($DATASET_SIZE * l($DATASET_SIZE))" | bc -l)
    floor_result=$(echo "$result / 1" | bc)
    timeout 3600 /usr/bin/time -v ${SNJ} -data true_aligned -seed 2 -n_i $floor_result -n_s $floor_log_n -n_o 3
    # cd ~/placement/test
    # python snj_postprocessing.py --dataset_size $n
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/SNJ/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt   
fi

if [ $run_type == "fnj" ]; then
    # Run fnj
    echo "Generating Newick tree file using fnj..."
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $FASTDIST -I fasta -O xml -D JC $MSA_FASTA_FILE > true_aligned.xml 
    timeout 1d /usr/bin/time -f "Time:\t%e s" $VMPEAK $FNJ -I xml -O newick true_aligned.xml  > ${RESULT_DIR}/fnj_msa.nwk 
    rm true_aligned.xml
fi