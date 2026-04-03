# examples/CMakeLists.txt
#
# Standalone build for the DIPPER library API examples.
# Build this separately from the main DIPPER project:
#
#   cd examples
#   mkdir build && cd build
#   cmake -DUSE_CUDA=ON ..
#   make -j$(nproc)
#
# The DIPPER source tree is referenced via DIPPER_SRC_DIR (default: ../src).
# Override it if your layout differs:
#   cmake -DUSE_CUDA=ON -DDIPPER_SRC_DIR=/absolute/path/to/dipper/src ..

cmake_minimum_required(VERSION 3.10)
project(dipper_examples LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-w)

option(USE_CUDA "Build with CUDA" OFF)

# Path to the DIPPER src/ directory (default: one level up from here).
set(DIPPER_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src"
    CACHE PATH "Path to the DIPPER src/ directory")

include(CheckLanguage)

# ---- CUDA setup --------------------------------------------------------------
set(BUILD_GPU_TARGET OFF)
if (USE_CUDA)
    check_language(CUDA)
    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        message(STATUS "CUDA compiler: ${CMAKE_CUDA_COMPILER}")
        set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86)
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -g -O3")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -O3 -rdc=true --extended-lambda -Xcompiler -std=c++17")
        set(BUILD_GPU_TARGET ON)
    else()
        message(FATAL_ERROR "USE_CUDA=ON but no CUDA compiler found.")
    endif()
else()
    message(FATAL_ERROR "Only CUDA is currently supported for the API examples. "
                        "Pass -DUSE_CUDA=ON to cmake.")
endif()

# ---- Dependencies ------------------------------------------------------------
add_compile_definitions(TBB_PREVIEW_GLOBAL_CONTROL=1)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB  REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})
find_package(TBB   REQUIRED CONFIG)

# ---- DIPPER source files shared across examples ------------------------------
set(DIPPER_SRCS
    ${DIPPER_SRC_DIR}/placement.cu
    ${DIPPER_SRC_DIR}/mash.cu
    ${DIPPER_SRC_DIR}/twoBitCompressor.cpp
    ${DIPPER_SRC_DIR}/fourBitCompressor.cpp
    ${DIPPER_SRC_DIR}/matrix_reader.cu
    ${DIPPER_SRC_DIR}/MSA.cu
    ${DIPPER_SRC_DIR}/tree.cpp
    ${DIPPER_SRC_DIR}/neighborJoining.cu
    ${DIPPER_SRC_DIR}/placement_close_k.cu
    ${DIPPER_SRC_DIR}/divide_and_conquer/mash.cu
    ${DIPPER_SRC_DIR}/divide_and_conquer/mash.cpp
    ${DIPPER_SRC_DIR}/divide_and_conquer/msa.cu
    ${DIPPER_SRC_DIR}/divide_and_conquer/placement_close_k.cpp
    ${DIPPER_SRC_DIR}/divide_and_conquer/placement_close_k.cu
)

# Helper macro: add a DIPPER example executable.
macro(add_dipper_example TARGET_NAME MAIN_SRC)
    add_executable(${TARGET_NAME} ${MAIN_SRC} ${DIPPER_SRCS})
    set_target_properties(${TARGET_NAME} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CXX_STANDARD              17
        CXX_STANDARD_REQUIRED     ON)
    target_include_directories(${TARGET_NAME} PRIVATE ${DIPPER_SRC_DIR})
    target_link_libraries(${TARGET_NAME} PRIVATE
        ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
    find_library(STDCXXFS NAMES stdc++fs)
    if (STDCXXFS)
        target_link_libraries(${TARGET_NAME} PRIVATE ${STDCXXFS})
    endif()
endmacro()

# ---- Example targets ---------------------------------------------------------
add_dipper_example(dipper_example_msa msa_placement_example.cpp)
