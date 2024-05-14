#include <cuda_runtime.h>
#include <iostream>

int main() {
    // Get device properties
    cudaDeviceProp prop;
    int device;

    // Get the default device ID
    cudaGetDevice(&device);

    // Get the properties of the device
    cudaGetDeviceProperties(&prop, device);

    std::cout << "Device Number: " << device << std::endl;
    std::cout << "Device Name: " << prop.name << std::endl;
    std::cout << "Maximum threads per block: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Maximum blocks per grid: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Maximum block dimensions: ("
              << prop.maxThreadsDim[0] << ", "
              << prop.maxThreadsDim[1] << ", "
              << prop.maxThreadsDim[2] << ")" << std::endl;
    std::cout << "Maximum grid dimensions: ("
              << prop.maxGridSize[0] << ", "
              << prop.maxGridSize[1] << ", "
              << prop.maxGridSize[2] << ")" << std::endl;
    std::cout << "Shared memory available per block: " << prop.sharedMemPerBlock << " bytes" << std::endl;

    return 0;
}