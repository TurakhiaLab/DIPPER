cmake_minimum_required(VERSION 3.10)

project(dipper LANGUAGES CXX)

option(USE_CUDA "Build with CUDA" OFF)
option(USE_HIP "Build with HIP" OFF)

if(NOT USE_CUDA AND NOT USE_HIP)
    message(FATAL_ERROR "Either USE_CUDA or USE_HIP must be enabled. e.g. cmake -DUSE_CUDA=ON ..")
endif()

include(CheckLanguage)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# ------------------ Versioning ------------------
set(dipper_VERSION "?.?.?")
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(dipper_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
configure_file(src/version.hpp.in ${CMAKE_SOURCE_DIR}/src/version.hpp)

# For CUDA
if(USE_CUDA)
        message(STATUS "Build with CUDA")
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
                enable_language(CUDA)
                message(STATUS "CUDA compiler detected at: ${CMAKE_CUDA_COMPILER}")
                TRY_RUN(
                        FIND_GPU_RESULT
                        FIND_GPU_COMPILE_RESULT
                        ${CMAKE_CURRENT_BINARY_DIR}/
                        ${CMAKE_CURRENT_BINARY_DIR}/../src/findgpu.cu
                        RUN_OUTPUT_VARIABLE GPU_ARCH_OUTPUT
                )
                if("${FIND_GPU_COMPILE_RESULT}" AND ("${FIND_GPU_RESULT}" EQUAL 0))
                        message(STATUS "CUDA-capable GPU found. GPU architecture: ${GPU_ARCH_OUTPUT}")
                        if(${CMAKE_VERSION} VERSION_LESS "3.18")
                                set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_${GPU_ARCH_OUTPUT}")
                        else()
                                set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCH_OUTPUT})
                        endif()
                else()
                        message(STATUS "No CUDA-capable GPU found. Generating binary to support multiple GPU generations.")
                        set(CMAKE_CUDA_ARCHITECTURES 75 80 86 90)
                endif()
        else()
                message(FATAL_ERROR "No CUDA found")
        endif()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O4 -Wall -Wextra")
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_STANDARD_REQUIRED ON)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -O4 -Wall -Wextra -rdc=true --extended-lambda -Xcompiler -std=c++17")

        set (SRC_FILES  src/tree_generation.cu
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
                        src/divide_and_conquer/placement_close_k.cu
                )

elseif(USE_HIP)
        message(STATUS "Build with HIP")
        find_package(HIP QUIET)
        if(HIP_FOUND)
                enable_language(HIP)
	        message(STATUS "HIP found: ${HIP_VERSION}")
		add_definitions(-D__HIP_PLATFORM_AMD__)
		set(CMAKE_CXX_COMPILER "/opt/rocm-${HIP_COMPILE_VERSION}/bin/hipcc")
		include_directories("/opt/rocm-${HIP_COMPILE_VERSION}/include/")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated")
		

                set (SRC_FILES  src/hip/tree_generation.cu.hip
                                src/hip/placement.cu.hip
                                src/hip/mash.cu.hip
                                src/twoBitCompressor.cpp
                                src/fourBitCompressor.cpp
                                src/hip/matrix_reader.cu.hip
                                src/hip/MSA.cu.hip
                                src/tree.cpp
                                src/hip/neighborJoining.cu.hip
                                src/hip/placement_close_k.cu.hip
                                src/hip/divide_and_conquer/mash.cu.hip
                                src/divide_and_conquer/mash.cpp
                                src/hip/divide_and_conquer/msa.cu.hip
                                src/divide_and_conquer/placement_close_k.cpp
                                src/hip/divide_and_conquer/placement_close_k.cu.hip
                )
        else()
	    message(FATAL_ERROR "No HIP found")
	endif()
else()
        message(FATAL_ERROR "Either USE_CUDA or USE_HIP must be enabled. e.g. cmake -DUSE_CUDA=ON ..")
endif()
# add_executable (dipper_cpu
#         src/tree_generation.cu
#         src/placement.cu
#         cpu/mash.cu
#         src/twoBitCompressor.cpp
#         src/fourBitCompressor.cpp
#         cpu/matrix_reader.cu
#         src/MSA.cu
#         src/tree.cpp
#         src/neighborJoining.cu
#         cpu/placement_close_k.cu
#         src/divide_and_conquer/mash.cu
#         src/divide_and_conquer/mash.cpp
#         src/divide_and_conquer/msa.cu
#         src/divide_and_conquer/placement_close_k.cpp
#         src/divide_and_conquer/placement_close_k.cu)

add_compile_options(-w)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})
find_package(TBB REQUIRED CONFIG)

add_executable(dipper ${SRC_FILES})

set_target_properties(dipper PROPERTIES CUDA_SEPARABLE_COMPILATION ON CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON) # set_target_properties(dipper_cpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
find_library(STDCXXFS NAMES stdc++fs)
if (STDCXXFS)
  target_link_libraries(dipper PRIVATE ${STDCXXFS}) #   target_link_libraries(dipper_cpu PRIVATE ${STDCXXFS})
endif()
target_link_libraries(dipper PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb) # target_link_libraries(dipper_cpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)

install(TARGETS dipper RUNTIME DESTINATION bin)
