# phylo-accel

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```

## Run Instructions
```
./test-llk ../test_data/gene1_msa.aln ../test_data/gene1.best.treefile 
./test-llk-gpu ../test_data/gene1_msa.aln ../test_data/gene1.best.treefile
```

