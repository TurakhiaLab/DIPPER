# Note of usage on Linux server:
# Copy "phylo-accel" to the home directory
# Copy the "placement folder" under "phylo-accel" to the home directory
# Copy "NeighborJoining-on-GPU" to the home directory
# And leave the "phylo-accel" under home directory as well
# And clone every baseline to the home directory
# And clone "phastsim" and "MAPLE" to the home directory
for n in 1000
do
    echo "-------------------Calculating $n-------------------------------------------"

    # Creating datasets

    # Creating dataset directory
    # mkdir /data/zec022/phastsim_datasets/dataset_$n 

    # Generating tree structure and branch lengths
    # To change average branch lengths, please go to line 11 in ~/phylo-accel/test/test.py
    # It is set to 3 sites per branch on average now
    # cd ~/phylo-accel/build/ 
    # ./gen_newick.out $n > ../test/rawtree.phy 
    # cd ~/phylo-accel/test
    # rm -r phastsim_dataset
    # python test.py

    # Generating sars-cov-2-like sequences
    # cd ~/phastSim
    # pip install -e .
    # python bin/phastSim --outpath /data/zec022/phastsim_datasets/dataset_$n/ --seed 8 --createFasta --createInfo --createNewick --treeFile ~/phylo-accel/test/modified_tree.phy --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon --reference phastSim/example/MN908947.3.fasta


    # Testing

    # Notes: 
    # Phastsim would store the aligned sequences (also the raw sequences) in /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta
    # This is a single file, however MashTree wants every sequence in separate files. Below is the splitting process
    # cd ~/phylo-accel/test
    # mkdir phastsim_dataset
    # mkdir /data/zec022/phastsim_datasets/dataset_$n/seqs
    # python test2.py --dataset_size $n
    # cd /data/zec022/phastsim_datasets/dataset_$n/seqs

    # Running MashTree by providing (raw seqs)
    # This can generate a distance matrix, but it is extremely slow, so not recommended for n>=2000
    # mkdir /data/zec022/mashtree/dataset_$n
    # mashtree --numcpus 12 *.fa --sketch-size 1000 --kmerlength 25 --outmatrix /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy --outtree /data/zec022/mashtree/dataset_$n/tree.nwk
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/mashtree/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    # Generating distance matrix by a parallelized program
    # time ~/phylo-accel/build/test-mash-gpu -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta > /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy
    

    # Running GPU accelerated version of conventional NJ by providing (distance matrix)
    # mkdir /data/zec022/NJ_GPU/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/NeighborJoining-on-GPU/build/test-nj -i /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy > /data/zec022/NJ_GPU/dataset_$n/tree.nwk
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/NJ_GPU/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    # Running RapidNJ by providing (distance matrix) or (MSA)
    # mkdir /data/zec022/rapidNJ/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/rapidNJ/bin/rapidnj /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -i pd -o t -x /data/zec022/rapidNJ/dataset_$n/tree.nwk
    # timeout 3600 /usr/bin/time -v ~/rapidNJ/bin/rapidnj /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i fa -o t -x /data/zec022/rapidNJ/dataset_$n/tree.nwk -c 32
    # python remove_extra_characters_in_rapidnj.py --dataset_size $n
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/rapidNJ/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    # Running Placement algorithm by providing (distance matrix) or (MSA) or (Raw seqs)
    # mkdir /data/zec022/placement/dataset_$n
    # cd /data/zec022/placement/dataset_$n
    # timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -i d -o t > /data/zec022/placement/dataset_$n/tree.nwk
    # timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i m -o t -t 2 > /data/zec022/placement/dataset_$n/tree.nwk
    # timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i r -o t > /data/zec022/placement/dataset_$n/tree.nwk
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/placement/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    # Running FastME by providing (distance matrix) or (MSA)
    # mkdir /data/zec022/fastme/dataset_$n
    # timeout 7200 /usr/bin/time -v ~/FastME/src/fastme -i /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -o /data/zec022/fastme/dataset_$n/tree.nwk
    # timeout 7200 /usr/bin/time -v ~/FastME/src/fastme -i /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -d K -o /data/zec022/fastme/dataset_$n/tree.nwk
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/fastme/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    # Running decenttree
    # mkdir /data/zec022/decenttree/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/decenttree/build/decenttree -in /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -out /data/zec022/decenttree/dataset_$n/tree.nwk -t NJ-V -nt 32
    # python ~/phylo-accel/test/test3.py
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/decenttree/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 

    mkdir /data/zec022/SNJ/dataset_$n
    python create_npy.py --dataset_size $n
    log_n=$(echo "l($n)" | bc -l)
    floor_log_n=$(echo "$log_n / 1" | bc)
    result=$(echo "scale=10; sqrt($n * l($n))" | bc -l)
    floor_result=$(echo "$result / 1" | bc)
    cd ~/SNJ
    timeout 3600 /usr/bin/time -v snj -data msa -seed 2 -n_i $floor_result -n_s $floor_log_n -n_o 3
    cd ~/placement/test
    python snj_postprocessing.py --dataset_size $n
    cd ~/MAPLE
    ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/SNJ/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    cat ~/temp.txt_RFdistances.txt 


done
