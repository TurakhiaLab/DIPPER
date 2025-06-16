cmake_minimum_required(VERSION 3.10)

project(msa LANGUAGES CXX CUDA)

# For CUDA
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
find_package(CUDA 10.0 REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_86 -rdc=true --extended-lambda")

find_package(Boost COMPONENTS program_options REQUIRED)

find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

# include(${TBB_DIR}/cmake/TBBBuild.cmake)
# tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
# find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

add_executable (dipper
                src/tree_generation.cu
                src/placement.cu
                src/mash.cu
                src/twoBitCompressor.cpp
                src/fourBitCompressor.cpp
                src/matrix_reader.cu
                src/MSA.cu
		src/tree.cpp
                src/neighborJoining.cu
                src/placement_close_k.cu
                src/divide_and_conquer/mash.cu
                src/divide_and_conquer/mash.cpp
                src/divide_and_conquer/msa.cu
                src/divide_and_conquer/placement_close_k.cpp
                src/divide_and_conquer/placement_close_k.cu)

find_package(TBB REQUIRED CONFIG)
set_target_properties(dipper PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(dipper PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)


