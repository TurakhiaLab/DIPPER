for n in 200
do
    echo "-------------------Calculating $n-------------------------------------------"
    cd ~/placement/test/apples-ii_testing
    ./convert_place < /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy > /data/zec022/phastsim_datasets/dataset_$n/apples_matrix.mat
    ./convert_backbone < /data/zec022/phastsim_datasets/dataset_$n/distance_matrix.phy > /data/zec022/phastsim_datasets/dataset_$n/apples_backbone.mat
    mkdir /data/zec022/apples-ii
    mkdir /data/zec022/apples-ii/dataset_$n
    ~/rapidNJ/bin/rapidnj /data/zec022/phastsim_datasets/dataset_$n/apples_backbone.mat -i pd -o t -x /data/zec022/apples-ii/dataset_$n/backbone.nwk
    cd ~/apples
    python run_apples.py -d /data/zec022/phastsim_datasets/dataset_$n/apples_matrix.mat -t /data/zec022/apples-ii/dataset_$n/backbone.nwk

    # python remove_bracket.py < output.nwk > processed.nwk
    # cd ~/MAPLE
    # ~/pypy3.10-v7.3.16-linux64/bin/pypy3 MAPLEv0.3.6.py --inputTree /data/zec022/phastsim_datasets/dataset_$n/sars-cov-2_simulation_output.tree  --inputRFtrees ~/placement/test/apples-ii_testing/processed.nwk --output ~/temp.txt --overwrite
    # cat ~/temp.txt_RFdistances.txt 


done
