#include <cuda_runtime_api.h>
#include <cusparseLt.h>       // cusparseLt header
#include <cstdio>             
#include <iostream>
#include <cmath>
#include "generateMatrix.hpp"  

int main(void) {
    // generate matrices for computation
    __half hA[N * N];                                                          
    __half hB[N * N];                                                          
    __half hC[N * N] = {};
    generate(hA, hB);
    float alpha = 1.0f;
    float beta  = 0.0f;
    const int runs = 100;
    //--------------------------------------------------------------------------
    // Device memory management
    //for (int iterations = 0; iterations < runs; iterations++) {
    __half *dA, *dB, *dC, *dD, *dA_compressed;
    int    *d_valid;
    cudaMalloc((void**) &dA, N*N*sizeof(__half));                           // allocate memory on the GPU (device) for the matrices
    cudaMalloc((void**) &dB, N*N*sizeof(__half));
    cudaMalloc((void**) &dC, N*N*sizeof(__half));
    cudaMalloc((void**) &d_valid, sizeof(int));
    dD = dC;

    cudaMemcpy(dA, hA, N*N*sizeof(__half), cudaMemcpyHostToDevice);           // copy host generated matrix to GPU (avoid page faulting)
    cudaMemcpy(dB, hB, N*N*sizeof(__half), cudaMemcpyHostToDevice);
    cudaMemcpy(dC, hC, N*N*sizeof(__half), cudaMemcpyHostToDevice);
    //--------------------------------------------------------------------------
    cusparseLtHandle_t             handle;                             // create handle and descriptors for the operation
    cusparseLtMatDescriptor_t      matA, matB, matC;
    cusparseLtMatmulDescriptor_t   matmul;
    cusparseLtMatmulAlgSelection_t alg_sel;
    cusparseLtMatmulPlan_t         plan;
    cudaStream_t                   stream = nullptr;
    cusparseLtInit(&handle);
    // matrix descriptor initialization
    cusparseLtStructuredDescriptorInit(                         
                                &handle, &matA, N,
                                N, N, 16,
                                CUDA_R_16F, CUSPARSE_ORDER_ROW,
                                CUSPARSELT_SPARSITY_50_PERCENT);
    cusparseLtDenseDescriptorInit(                              
                                &handle, &matB, N,
                                N, N, 16,
                                CUDA_R_16F, CUSPARSE_ORDER_ROW);
    cusparseLtDenseDescriptorInit(                              
                                &handle, &matC, N,
                                N, N, 16,
                                CUDA_R_16F, CUSPARSE_ORDER_ROW);
    // matmul, algorithm selection, and plan initialization
    cusparseLtMatmulDescriptorInit(                            
                                &handle, &matmul, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &matA, &matB, &matC, &matC,
                                CUSPARSE_COMPUTE_16F);
    cusparseLtMatmulAlgSelectionInit(                          
                                &handle, &alg_sel, &matmul,
                                CUSPARSELT_MATMUL_ALG_DEFAULT);
    cusparseLtMatmulPlanInit(   
                                &handle, &plan, &matmul, &alg_sel);


    // Compress the A matrix
    size_t compressed_size, compressed_buffer_size;                            
    void*  dA_compressedBuffer;
    cusparseLtSpMMACompressedSize(
                                &handle, &plan,
                                &compressed_size,
                                &compressed_buffer_size);
    cudaMalloc((void**) &dA_compressed, compressed_size);
    cudaMalloc((void**) &dA_compressedBuffer, compressed_buffer_size);

    cusparseLtSpMMACompress(&handle, &plan, dA, dA_compressed, dA_compressedBuffer,stream);

    // see structure of compressed matrix
    __half hA_compressed[compressed_size];
    __half hA_compressed_buffer[compressed_buffer_size];
    cudaMemcpy(hA_compressed, dA_compressed, compressed_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(hA_compressed_buffer, dA_compressedBuffer, compressed_buffer_size, cudaMemcpyDeviceToHost);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Search the best kernel
    int           num_streams = 0;                                              // search API for the best kernel for the matmul operation
    cudaStream_t* streams     = nullptr;
    cusparseLtMatmulSearch(
                                &handle, &plan, &alpha,
                                dA_compressed, dB, &beta,
                                dC, dD, nullptr,
                                streams, num_streams);

    size_t workspace_size;
    cusparseLtMatmulPlanInit(&handle, &plan, &matmul, &alg_sel); // redundant?

    cusparseLtMatmulGetWorkspace(
                                &handle, &plan,                
                                &workspace_size);
    void* d_workspace;
    cudaMalloc((void**) &d_workspace, workspace_size);
    //for (int iterations = 0; iterations < runs; iterations++) {
        // Perform the matrix multiplication
        cusparseLtMatmul(
                                    &handle, &plan, &alpha, 
                                    dA_compressed, dB,    // do the actual multiplication operation finally, in the allocated workspace
                                    &beta, dC, dD, d_workspace, streams,
                                    num_streams);
    //}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // destroy plan and handle                                          
    cusparseLtMatDescriptorDestroy(&matA);                 // deallocate no longer necessary handles
    cusparseLtMatDescriptorDestroy(&matB);
    cusparseLtMatDescriptorDestroy(&matC);
    cusparseLtMatmulPlanDestroy(&plan);
    cusparseLtDestroy(&handle);
    
    // //--------------------------------------------------------------------------
    // // device result check
    // // matrix A has been pruned
    // cudaMemcpy(hA, dA, N*N*sizeof(__half), cudaMemcpyDeviceToHost);            // check that computation worked.
    // cudaMemcpy(hC, dC, N*N*sizeof(__half), cudaMemcpyDeviceToHost);

    // // host computation
    // float hC_result[N * N];
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         float sum  = 0.0f;
    //         for (int k1 = 0; k1 < N; k1++) {
    //             auto posA = i * N + k1;
    //             auto posB = k1 * N + j;
    //             sum      += static_cast<float>(hA[posA]) *  // [i][k]
    //                         static_cast<float>(hB[posB]);   // [k][j]
    //         }
    //         auto posC       = i * N + j;
    //         hC_result[posC] = sum;  // [i][j]
    //     }
    // }
    // // host-device comparison
    // int correct = 1;
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         auto pos          = i * N + j;
    //         auto device_value = static_cast<float>(hC[pos]);
    //         auto host_value   = hC_result[pos];
    //         if (std::abs(device_value - host_value) > 2) {
    //             // direct floating point comparison is not reliable
    //             std::printf("(%d, %d):\t%f vs. %f\n",
    //                         i, j, host_value, device_value);
    //             correct = 0;
    //             break;
    //         }
    //     }
    // }
    // if (correct)
    //     std::printf("matmul_example test PASSED\n");
    // else
    //     std::printf("matmul_example test FAILED: wrong result\n");

    // device memory deallocation
    cudaFree(dA_compressed);
    cudaFree(dA);
    cudaFree(dB);
    cudaFree(dC);
    cudaFree(d_valid);
    cudaFree(d_workspace);
    cudaFree(dA_compressedBuffer);
    //}
    return EXIT_SUCCESS;
}
