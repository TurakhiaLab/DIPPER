#include <cuda_runtime_api.h> // cudaMalloc, cudaMemcpy, etc.
#include <cusparseLt.h>       // cusparseLt header
#include <cstdio>             // printf
#include <cstdlib>            // std::rand
#include "cudaErrorHandling.h" // error status handling


constexpr int EXIT_UNSUPPORTED = 2;

int main(void) {
    // Host problem definition, row-major order
    // bigger sizes may require dynamic allocations
    constexpr int m            = 32;                                            // columns
    constexpr int n            = 32;                                            // rows
    constexpr int k            = 32;                                            // inner dimension
    auto          order        = CUSPARSE_ORDER_ROW;                            // specifies data is in row major format. other formats: column major
    auto          opA          = CUSPARSE_OPERATION_NON_TRANSPOSE;              // specify matrices are not transposed in the operation 
    auto          opB          = CUSPARSE_OPERATION_NON_TRANSPOSE;
    auto          type         = CUDA_R_16F;                                    // data type is float16
    auto          compute_type = CUSPARSE_COMPUTE_16F;                          // compute type for cusparse is also float16. What is difference between CUDA_R_16F and CUSPARSE_COMPUTE_16F?

    bool     is_rowmajor    = (order == CUSPARSE_ORDER_ROW);                   // various information about the matrix operation to be done. These values referenced throughout by other functions. 
    bool     isA_transposed = (opA != CUSPARSE_OPERATION_NON_TRANSPOSE);
    bool     isB_transposed = (opB != CUSPARSE_OPERATION_NON_TRANSPOSE);
    auto     num_A_rows     = (isA_transposed) ? k : m;
    auto     num_A_cols     = (isA_transposed) ? m : k;
    auto     num_B_rows     = (isB_transposed) ? n : k;
    auto     num_B_cols     = (isB_transposed) ? k : n;
    auto     num_C_rows     = m;
    auto     num_C_cols     = n;
    unsigned alignment      = 16;                                              // bytes of each element of matrix in memory
    auto     lda            = (is_rowmajor) ? num_A_cols : num_A_rows;
    auto     ldb            = (is_rowmajor) ? num_B_cols : num_B_rows;
    auto     ldc            = (is_rowmajor) ? num_C_cols : num_C_rows;
    auto     A_height       = (is_rowmajor) ? num_A_rows : num_A_cols;
    auto     B_height       = (is_rowmajor) ? num_B_rows : num_B_cols;
    auto     C_height       = (is_rowmajor) ? num_C_rows : num_C_cols;
    auto     A_size         = A_height * lda * sizeof(__half);
    auto     B_size         = B_height * ldb * sizeof(__half);
    auto     C_size         = C_height * ldc * sizeof(__half);
    __half hA[m * k];                                                          
    __half hB[k * n];                                                          // half is just float16
    __half hC[m * n] = {};
    for (int i = 0; i < m * k; i++)                                            // randomly generate two matrices
        hA[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
    for (int i = 0; i < k * n; i++)
        hB[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
    float alpha = 1.0f;
    float beta  = 0.0f;
    //--------------------------------------------------------------------------
    // Device memory management
    __half *dA, *dB, *dC, *dD, *dA_compressed;
    int    *d_valid;
    CHECK_CUDA( cudaMalloc((void**) &dA, A_size) )                             // allocate memory on the GPU (device) for the matrices
    CHECK_CUDA( cudaMalloc((void**) &dB, B_size) )
    CHECK_CUDA( cudaMalloc((void**) &dC, C_size) )
    CHECK_CUDA( cudaMalloc((void**) &d_valid, sizeof(int)) )
    dD = dC;

    CHECK_CUDA( cudaMemcpy(dA, hA, A_size, cudaMemcpyHostToDevice) )            // copy host generated matrix to GPU (avoid page faulting)
    CHECK_CUDA( cudaMemcpy(dB, hB, B_size, cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(dC, hC, C_size, cudaMemcpyHostToDevice) )
    //--------------------------------------------------------------------------
    cusparseLtHandle_t             handle;                                      // create handle and descriptors for the operation
    cusparseLtMatDescriptor_t      matA, matB, matC;
    cusparseLtMatmulDescriptor_t   matmul;
    cusparseLtMatmulAlgSelection_t alg_sel;
    cusparseLtMatmulPlan_t         plan;
    cudaStream_t                   stream = nullptr;
    CHECK_CUSPARSE( cusparseLtInit(&handle) )
    // matrix descriptor initialization
    CHECK_CUSPARSE( cusparseLtStructuredDescriptorInit(                         // create descriptor for matA using information about it (from above)
                                            &handle, &matA, num_A_rows,
                                            num_A_cols, lda, alignment,
                                            type, order,
                                            CUSPARSELT_SPARSITY_50_PERCENT) )
    CHECK_CUSPARSE( cusparseLtDenseDescriptorInit(                              // create descriptor for matB using information about it (from above)
                                            &handle, &matB, num_B_rows,
                                            num_B_cols, ldb, alignment,
                                            type, order) )
    CHECK_CUSPARSE( cusparseLtDenseDescriptorInit(                              // create descriptor for matC using information about it (from above)
                                            &handle, &matC, num_C_rows,
                                            num_C_cols, ldc, alignment,
                                            type, order) )
    // matmul, algorithm selection, and plan initialization
    CHECK_CUSPARSE( cusparseLtMatmulDescriptorInit(                            // create descriptor for the operation (matmul)
                                            &handle, &matmul, opA, opB,
                                            &matA, &matB, &matC, &matC,
                                            compute_type) )
    CHECK_CUSPARSE( cusparseLtMatmulAlgSelectionInit(                          // select multiplication algorithm (default). What are other algorithms and how does performance compare?
                                            &handle, &alg_sel, &matmul,
                                            CUSPARSELT_MATMUL_ALG_DEFAULT) )    // CUSPARSELT_MATMUL_ALG_DEFAULT is the only possible value according to documentation. Why make it a flag?
    CHECK_CUSPARSE( cusparseLtMatmulPlanInit(&handle, &plan, &matmul, &alg_sel))  // returns the size of the workspace and the compiled kernel needed

    //--------------------------------------------------------------------------
    // Prune the A matrix (in-place) and check the correctness
    CHECK_CUSPARSE( cusparseLtSpMMAPrune(&handle, &matmul, dA, dA,             // unnecessary if input matrix already satisfies 2:4 sparsity
                                         CUSPARSELT_PRUNE_SPMMA_TILE, stream) )
    CHECK_CUSPARSE( cusparseLtSpMMAPruneCheck(&handle, &matmul, dA,
                                              d_valid, stream) )
    int is_valid;
    CHECK_CUDA( cudaMemcpyAsync(&is_valid, d_valid, sizeof(int),
                                cudaMemcpyDeviceToHost, stream) )
    CHECK_CUDA( cudaStreamSynchronize(stream) )
    if (is_valid != 0) {
        std::printf("!!!! The matrix has been pruned in a wrong way. "
                    "cusparseLtMatmul will not provide correct results\n");
        return EXIT_FAILURE;
    }
    //--------------------------------------------------------------------------
    // Compress the A matrix
    size_t compressed_size, compressed_buffer_size;                             // compress matrix, following plan
    void*  dA_compressedBuffer;
    CHECK_CUSPARSE( cusparseLtSpMMACompressedSize(&handle, &plan,
                                                  &compressed_size,
                                                  &compressed_buffer_size) )
    CHECK_CUDA( cudaMalloc((void**) &dA_compressed, compressed_size) )
    CHECK_CUDA( cudaMalloc((void**) &dA_compressedBuffer,
                           compressed_buffer_size) )

    CHECK_CUSPARSE( cusparseLtSpMMACompress(&handle, &plan, dA, dA_compressed,
                                            dA_compressedBuffer,stream) )
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Search the best kernel
    int           num_streams = 0;                                              // search API for the best kernel for the matmul operation
    cudaStream_t* streams     = nullptr;
    CHECK_CUSPARSE( cusparseLtMatmulSearch(&handle, &plan, &alpha,
                                           dA_compressed, dB, &beta,
                                           dC, dD, nullptr,
                                           streams, num_streams) )
    // otherwise, it is possible to set it directly:
    //int alg = 0;
    //CHECK_CUSPARSE( cusparseLtMatmulAlgSetAttribute(
    //                                        &handle, &alg_sel,
    //                                        CUSPARSELT_MATMUL_ALG_CONFIG_ID,
    //                                        &alg, sizeof(alg)))
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    size_t workspace_size;
    CHECK_CUSPARSE( cusparseLtMatmulPlanInit(&handle, &plan, &matmul, &alg_sel)) // redundant

    CHECK_CUSPARSE( cusparseLtMatmulGetWorkspace(&handle, &plan,                // fetch plan, workspace size as defined by PlanInit
                                                 &workspace_size))
    void* d_workspace;
    CHECK_CUDA( cudaMalloc((void**) &d_workspace, workspace_size) )
    // Perform the matrix multiplication
    CHECK_CUSPARSE( cusparseLtMatmul(&handle, &plan, &alpha, dA_compressed, dB,    // do the actual multiplication operation finaly, in the allocated workspace
                                     &beta, dC, dD, d_workspace, streams,
                                     num_streams) )
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // destroy plan and handle                                          
    CHECK_CUSPARSE( cusparseLtMatDescriptorDestroy(&matA) )                     // deallocate no longer necessary handles
    CHECK_CUSPARSE( cusparseLtMatDescriptorDestroy(&matB) )
    CHECK_CUSPARSE( cusparseLtMatDescriptorDestroy(&matC) )
    CHECK_CUSPARSE( cusparseLtMatmulPlanDestroy(&plan) )
    CHECK_CUSPARSE( cusparseLtDestroy(&handle) )
    //--------------------------------------------------------------------------
    // device result check
    // matrix A has been pruned
    CHECK_CUDA( cudaMemcpy(hA, dA, A_size, cudaMemcpyDeviceToHost) )            // use pruned version of matrix to check that computation worked.
    CHECK_CUDA( cudaMemcpy(hC, dC, C_size, cudaMemcpyDeviceToHost) )

    bool A_std_layout = (is_rowmajor != isA_transposed);
    bool B_std_layout = (is_rowmajor != isB_transposed);
    // host computation
    float hC_result[m * n];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            float sum  = 0.0f;
            for (int k1 = 0; k1 < k; k1++) {
                auto posA = (A_std_layout) ? i * lda + k1 : i + k1 * lda;
                auto posB = (B_std_layout) ? k1 * ldb + j : k1 + j * ldb;
                sum      += static_cast<float>(hA[posA]) *  // [i][k]
                            static_cast<float>(hB[posB]);   // [k][j]
            }
            auto posC       = (is_rowmajor) ? i * ldc + j : i + j * ldc;
            hC_result[posC] = sum;  // [i][j]
        }
    }
    // host-device comparison
    int correct = 1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            auto pos          = (is_rowmajor) ? i * ldc + j : i + j * ldc;
            auto device_value = static_cast<float>(hC[pos]);
            auto host_value   = hC_result[pos];
            if (device_value != host_value) {
                // direct floating point comparison is not reliable
                std::printf("(%d, %d):\t%f vs. %f\n",
                            i, j, host_value, device_value);
                correct = 0;
                break;
            }
        }
    }
    if (correct)
        std::printf("matmul_example test PASSED\n");
    else
        std::printf("matmul_example test FAILED: wrong result\n");
    //--------------------------------------------------------------------------
    // device memory deallocation
    CHECK_CUDA( cudaFree(dA_compressed) )
    CHECK_CUDA( cudaFree(dA) )
    CHECK_CUDA( cudaFree(dB) )
    CHECK_CUDA( cudaFree(dC) )
    CHECK_CUDA( cudaFree(d_valid) )
    CHECK_CUDA( cudaFree(d_workspace) )
    CHECK_CUDA( cudaFree(dA_compressedBuffer) )
    return EXIT_SUCCESS;
}
