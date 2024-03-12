# ECE 284 FINAL PROJECT
# FAST DISTANCE ESTIMATION USING MASH

# SERIAL IMPLEMENTATION 
## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make test-mash-gpu-serial
```

## Run Instructions to reproduce key results
```
../run_commands_gpu_serial.sh
```

# PARALLEL IMPLEMENTATION
## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make test-mash-gpu
```

## Run Instructions to reproduce key results
```
../run_commands_gpu_parallel.sh
```

## Instructions to extract key results from the runs

#  Mash Distance Computation Kernel (mashDistConstruction)

   # Result 1: Speedup vs numSequences

    Fixed Parameters:

        sketchSize = 1000
        kmerSize = 15
        numBlocks = 128

    Baseline performance numbers for serial GPU implementation will be present in: <git directory>/run_logs/serial_sars<numSequences>.txt files

    Performance numbers for parallel GPU implementation will be present in: <git directory>/run_logs/parallel_sars<numSequences>.txt files

    Please refer to run_commands_gpu_parallel.sh and run_commands_gpu_serial.sh to see which results are present under which file names.

    Look for profiled performance numbers for mashDistConstruction kernel in the above files to find speedup corresponding to various values of numSequences

   # Result 2: Speedup vs numBlocks

    Fixed Parameters:

        numSequences = 2000
        sketchSize = 1000
        kmerSize = 15

    Baseline performance numbers for serial GPU implementation will be present in: <git directory>/run_logs/serial_sars2000.txt 

    Performance numbers for parallel GPU implementation will be present in: <git directory>/run_logs/parallel_sars2000_b<blockSize>.txt files

    Look for profiled performance numbers for mashDistConstruction kernel in the above files to find speedup corresponding to various values of numBlocks


   # Result 3: Speedup vs sketchSize

       Fixed Parameters:

        numSequences = 2000
        numBlocks = 128
        kmerSize = 15
        
    Baseline performance numbers for serial GPU implementation will be present in: <git directory>/run_logs/serial_sars2000_s<sketchSize>.txt 

    Performance numbers for parallel GPU implementation will be present in: <git directory>/run_logs/parallel_sars2000_s<sketchSize>.txt files

    Look for profiled performance numbers for mashDistConstruction kernel in the above files to find speedup corresponding to various values of numBlocks


# Sketch Computation Kernel (sketchConstruction)

   # Result 1: Speedup vs numBlocks

    Fixed Parameters:

        numSequences = 2000
        sketchSize = 1000
        kmerSize = 15

    Baseline performance numbers for serial GPU implementation will be present in: <git directory>/run_logs/serial_sars2000.txt 

    Performance numbers for parallel GPU implementation will be present in: <git directory>/run_logs/parallel_sars2000_b<blockSize>.txt files

    Look for profiled performance numbers for sketchConstruction kernel in the above files to find speedup corresponding to various values of numBlocks

   # Result 2: Speedup vs numSequences

    Fixed Parameters:

        sketchSize = 1000
        kmerSize = 15
        numBlocks = 128

    Baseline performance numbers for serial GPU implementation will be present in: <git directory>/run_logs/serial_sars<numSequences>.txt files

    Performance numbers for parallel GPU implementation will be present in: <git directory>/run_logs/parallel_sars<numSequences>.txt files

    Please refer to run_commands_gpu_parallel.sh and run_commands_gpu_serial.sh to see which results are present under which file names.

    Look for profiled performance numbers for sketchConstruction kernel in the above files to find speedup corresponding to various values of numSequences