
# DIPPER-TWILIGHT
results_dir=/data/swalia/dipper/alisim/results/iterative/dipper-twilight-10000
ref_tree=/data/swalia/dipper/alisim/trees
mkdir -p maple
touch maple/nrf

l="10000 20000 50000 100000 200000 500000 1000000"

for d in $l; do
    #ref_alignment=/data/swalia/dipper/alisim/datasets/${d}/len_10000_leaves_${d}.aligned.fa
    ref_seqs=/data/swalia/dipper/alisim/datasets/${d}/len_10000_leaves_${d}.unaligned.fa

    # /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh dipper -f ${ref_seqs} -s 1000 -t 2 -i r -o t > ${results_dir}/dipper_${d}.nwk
    /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh /home/swalia@AD.UCSD.EDU/work/TWILIGHT/build/twilight -t ${results_dir}/dipper_${d}.nwk -i ${ref_seqs} -o ${results_dir}/twilight_${d}.aln -C 32 -v  --length-deviation 1 -G 1 
    # python3 /home/swalia@AD.UCSD.EDU/tools/MAPLE/MAPLEv0.3.6.py --inputTree ${ref_tree}/nleaves${d}.treefile --inputRFtrees ${results_dir}/dipper_${d}.nwk --output maple/nrf --overwrite
    # echo "nrf iteration $itr $(awk 'NR==2 {print $2}' maple/nrf_RFdistances.txt)"

    # FastSP.jar -r ${ref_alignment} -e ${results_dir}/twilight_${d}.aln
done
# 
# itr=0
# /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh dipper -f ${ref_seqs} -s 1000 -t 2 -i r -o t > ${results_dir}/dipper_itr_${itr}.nwk
# /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh /home/swalia@AD.UCSD.EDU/work/TWILIGHT/build/twilight -t ${results_dir}/dipper_itr_${itr}.nwk -i ${ref_seqs} -o ${results_dir}/twilight_itr_${itr}.aln -C 32 -v  --length-deviation 1

# prev_itr=0
# for itr in {1..5}; do

#     /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh dipper -f ${results_dir}/twilight_itr_${prev_itr}.aln -s 1000 -t 2 -i m -o t > ${results_dir}/dipper_itr_${itr}.nwk
#     /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh /home/swalia@AD.UCSD.EDU/work/TWILIGHT/build/twilight -t ${results_dir}/dipper_itr_${itr}.nwk -i ${ref_seqs} -o ${results_dir}/twilight_itr_${itr}.aln -C 32 -v  --length-deviation 1
#     prev_itr=${itr}
# done

# for itr in {0..5}; do
#      python3 /home/swalia@AD.UCSD.EDU/tools/MAPLE/MAPLEv0.3.6.py --inputTree ${ref_tree} --inputRFtrees ${results_dir}/dipper_itr_${itr}.nwk --output maple/nrf --overwrite
#      echo "nrf iteration $itr $(awk 'NR==2 {print $2}' maple/nrf_RFdistances.txt)"
# done

# # Pasttree-TWILIGHT-Fasttree
# results_dir=/data/swalia/dipper/alisim/results/iterative/parttree-twilight-fasttree-10000
# ref_tree=/data/swalia/dipper/alisim/trees/nleaves10000_iterative.treefile
# ref_alignment=/data/swalia/dipper/alisim/datasets/iterative/10000/10000.aligned.fa
# ref_seqs=/data/swalia/dipper/alisim/datasets/iterative/10000/10000.unaligned.fa

# /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh mafft --retree 0 --treeout --parttree --reorder --thread 32 ${ref_seqs} > mafft.out
# rm mafft.out
# mv ${ref_seqs}.tree parttree.tree
# python3 /data/swalia/dipper/alisim/results/iterative/mafft2nwk.py parttree.tree ${ref_seqs} parttree_itr_0.nwk --parttree
# rm parttree.tree
# /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh /home/swalia@AD.UCSD.EDU/work/TWILIGHT/build/twilight -t ${results_dir}/parttree_itr_0.nwk -i ${ref_seqs} -o ${results_dir}/twilight_itr_0.aln -C 32 -v --no-filtering

# prev_itr=0
# for itr in {1..5}; do
#     /usr/bin/time -f "Time:\t%e s" ~/work/phylo-accel/scripts/vmPeak.sh fasttree -nt -nocat -noml twilight_itr_${prev_itr}.aln  > fasttree_itr_${itr}.nwk
#     /home/swalia@AD.UCSD.EDU/work/TWILIGHT/build/twilight -t fasttree_itr_${itr}.nwk -i ${ref_seqs} -o twilight_itr_${itr}.aln -C 32 --no-filtering
#     prev_itr=${itr}
# done

# python3 /home/swalia@AD.UCSD.EDU/tools/MAPLE/MAPLEv0.3.6.py --inputTree ${ref_tree} --inputRFtrees ${results_dir}/parttree_itr_0.nwk --output maple/nrf --overwrite
# echo "nrf iteration 0 $(awk 'NR==2 {print $2}' maple/nrf_RFdistances.txt)"

# for itr in {1..2}; do
#     python3 /home/swalia@AD.UCSD.EDU/tools/MAPLE/MAPLEv0.3.6.py --inputTree ${ref_tree} --inputRFtrees ${results_dir}/fasttree_itr_${itr}.nwk --output maple/nrf --overwrite
#     echo "nrf iteration $itr $(awk 'NR==2 {print $2}' maple/nrf_RFdistances.txt)"
# done
