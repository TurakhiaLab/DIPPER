# phylo-accel

## Make sure going to placement subfolder. Everything is there.

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make mash-placement
```

## Run Instructions
```
./mash-placement -f [input-file] -i [r or m or d] -o t > [output-file]
# -i r -> raw sequences
# -i m -> MSA
# -i d -> distance matrix
# for more options, try ./mash-placement -h
```

