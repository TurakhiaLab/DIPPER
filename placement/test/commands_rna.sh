for n in 100000
do
    echo "-------------------Calculating $n-------------------------------------------"
    # mkdir /data/zec022/phastsim_datasets/dataset_$n
    # cd ~/phylo-accel/build/
    # ./gen_newick.out $n > ../test/rawtree.phy
    # cd ~/phylo-accel/test
    # rm -r phastsim_dataset
    # python test.py
    # cd ~/phastSim
    # python bin/phastSim --outpath /data/zec022/phastsim_datasets/dataset_$n/ --seed 8 --createFasta --createInfo --createNewick --treeFile ~/phylo-accel/test/modified_tree.phy --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon --reference phastSim/example/MN908947.3.fasta

    # cd ~/phylo-accel/test
    # mkdir phastsim_dataset
    # mkdir /data/zec022/seqs
    # python test2.py
    # cd /data/zec022/seqs
    # mashtree --numcpus 12 *.fa --sketch-size 1000 --kmerlength 25 --outmatrix /data/zec022/distance_matrix.phy --outtree ~/phylo-accel/test/phastsim_dataset/nj_tree.nwk 
    # python ~/phylo-accel/test/test4.py
    # mkdir /data/zec022/phastsim_datasets/dataset_$n
    # time ~/phylo-accel/build/test-mash-gpu -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta > /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy
    # mkdir /data/zec022/NJ_GPU/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/NeighborJoining-on-GPU/build/test-nj -i /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy > /data/zec022/NJ_GPU/dataset_$n/tree.nwk
    # mkdir /data/zec022/rapidNJ/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/rapidNJ/bin/rapidnj /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -i pd -o t -x /data/zec022/rapidNJ/dataset_$n/tree.nwk
    # timeout 3600 /usr/bin/time -v ~/rapidNJ/bin/rapidnj /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i fa -o t -x /data/zec022/rapidNJ/dataset_$n/tree.nwk -c 32
    # mkdir /data/zec022/placement/dataset_$n
    # cd /data/zec022/placement/dataset_$n
    # timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -i d -o t 
    # > /data/zec022/placement/dataset_$n/tree.nwk
    # timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i m -o t 
    timeout 3600 /usr/bin/time -v  ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i m -o t -t 2 > /data/zec022/placement/dataset_$n/tree.nwk
    # rm /data/zec022/placement/dataset_$n/analysis_report*
    # timeout 3600 /usr/bin/time -v  nsys profile -o analysis_report --force-overwrite=true ~/placement/build/mash-placement -f /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.fasta -i m -o t > /data/zec022/placement/dataset_$n/tree.nwk
    # nsys stats analysis_report.nsys-rep
    # mkdir /data/zec022/fastme/dataset_$n
    # timeout 7200 /usr/bin/time -v ~/FastME/src/fastme -i /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -o /data/zec022/fastme/dataset_$n/tree.nwk
    # mkdir /data/zec022/decenttree/dataset_$n
    # timeout 3600 /usr/bin/time -v ~/decenttree/build/decenttree -in /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy -out /data/zec022/decenttree/dataset_$n/tree.nwk -t NJ-V -nt 32
    # python ~/phylo-accel/test/test3.py
    # timeout 3600 /usr/bin/time -v  ~/phylo-accel/build/placement < /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy > /data/zec022/placement/dataset_$n/tree.nwk

    cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/NJ_GPU/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 
    # sed -i "s/'//g" /data/zec022/rapidNJ/dataset_$n/tree.nwk
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/rapidNJ/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/placement/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/decenttree/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees /data/zec022/fastme/dataset_$n/tree.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 


done
