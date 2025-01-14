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

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)


add_executable(generateDis test/generateDistanceMatrix.cpp)
add_executable (msa
                src/msa.cpp
                src/tree.cpp
                src/align.cpp)

add_executable (test-align
                test/test_align.cpp
                src/tree.cpp
                src/align.cpp)

add_executable (test-msa
                test/test_msa.cpp
                src/tree.cpp
                src/align.cpp)


add_executable (test-mash
                src/murmurHash3.cpp
                src/createSketch.cpp
                src/mashDistance.cpp
                src/distanceMatrix.cpp
                src/twoBitCompressor.cpp
                src/tree.cpp
                test/test_mash.cpp)

add_executable (test-mash-gpu
                src/mash.cu
                src/neighbourJoining.cu
                src/twoBitCompressor.cpp
                test/test_mash_gpu.cpp)

set_target_properties(test-mash-gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

                

add_executable (test-nj
                src/mash.cu
                src/rapidNJ.cu
                src/twoBitCompressor.cpp
                test/test_nj.cpp)

add_executable (placement
                src/test_placement.cpp
                src/placement.cu)

set_target_properties(placement PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable (mash-placement
                src/test_mash_placement.cpp
                src/mash_placement.cu
                src/twoBitCompressor.cpp)

set_target_properties(mash-placement PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


target_link_libraries(msa PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
target_link_libraries(test-align PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
target_link_libraries(test-msa PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
target_link_libraries(test-mash PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
target_link_libraries(test-mash-gpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
target_link_libraries(test-nj PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
target_link_libraries(generateDis PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
target_link_libraries(placement PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
target_link_libraries(mash-placement PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})

