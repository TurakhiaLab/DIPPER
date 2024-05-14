#include <iostream>
#include <cub/cub.cuh>

__global__ void testKernel1() {
    // Specialize BlockRadixSort for a 1D block of 10 threads owning 2 integer items each
    typedef cub::BlockRadixSort<uint64_t, 10, 2> BlockRadixSort;

    __shared__ typename BlockRadixSort::TempStorage temp_storage;

    // Make sure each thread has unique keys for clear sorting
    uint64_t keys[] = {threadIdx.x*3, threadIdx.x*2 + 1};

    BlockRadixSort(temp_storage).Sort(keys);

    printf("Thread %d has key[0]: %ld, key[1]: %ld\n", threadIdx.x, keys[0], keys[1]);
}

int main() {
    testKernel1<<<1,10,1000>>>();

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
        return 1;
    }
    cudaDeviceSynchronize();
    return 0;
}
