Files used for benchmarking performance of cuSPARSELt for Dense-Sparse Matrix Multiplication

Files can be compiled on peregrine using:

```
nvcc -o <output_file> <input_file> -arch=native
```
Additionally, flag "--keep" can be specified to make sure intermediate compiler files are kept in local directory instead of /tmp

To generate benchmark reports:
```
nsys profile <application>
```
This creates a .nsys-rep file that can be opened with [NSight Systems][https://developer.nvidia.com/nsight-systems].

