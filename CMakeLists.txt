cmake_minimum_required(VERSION 3.10)

project(dipper LANGUAGES CXX)

option(USE_CUDA "Build with CUDA" OFF)
option(USE_HIP "Build with HIP" OFF)
option(USE_CPU "Build with CPU" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT USE_CPU AND NOT USE_CUDA AND NOT USE_HIP)
    message("Building failed. Set either USE_CPU/USE_CUDA/USE_HIP")
    exit()
endif()
if(NOT USE_CUDA AND NOT USE_HIP)
    message("Building CPU only version. To build with CUDA, set USE_CUDA to ON. To build with HIP, set USE_HIP to ON.")
endif()

add_compile_options(-w)

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

set(BUILD_GPU_TARGET OFF)

# For CUDA
if(USE_CUDA)
        message(STATUS "Build with CUDA")
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
                enable_language(CUDA)
                message(STATUS "CUDA compiler detected at: ${CMAKE_CUDA_COMPILER}")
                # Optionally allow the user to set CUDA architectures manually, otherwise set a default
                if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
                    # Set default architectures (compatible for most modern NVIDIA GPUs)
                    set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80 86)
                endif()
                message(STATUS "Using CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
                set(BUILD_GPU_TARGET ON)
        else()
                message(STATUS "No CUDA found")
        endif()


        if(BUILD_GPU_TARGET)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O4 -Wall -Wextra")
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
        endif()
# For ROCM
elseif(USE_HIP)
        message(STATUS "Build with HIP")
        find_package(HIP QUIET)
        if(HIP_FOUND)
                enable_language(HIP)
	        message(STATUS "HIP found: ${HIP_VERSION}")
		add_definitions(-D__HIP_PLATFORM_AMD__)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated")
		

                set(HIP_GPU_AVAILABLE OFF)
                find_program(ROCMINFO_EXECUTABLE rocminfo)
                if(ROCMINFO_EXECUTABLE)
                        execute_process(COMMAND ${ROCMINFO_EXECUTABLE} RESULT_VARIABLE ROCMINFO_RESULT OUTPUT_QUIET ERROR_QUIET)
                        if(ROCMINFO_RESULT EQUAL 0)
                                set(HIP_GPU_AVAILABLE ON)
                        endif()
                endif()
                if(NOT HIP_GPU_AVAILABLE)
                        find_program(ROCMSMI_EXECUTABLE rocm-smi)
                        if(ROCMSMI_EXECUTABLE)
                                execute_process(COMMAND ${ROCMSMI_EXECUTABLE} -i RESULT_VARIABLE ROCMSMI_RESULT OUTPUT_QUIET ERROR_QUIET)
                                if(ROCMSMI_RESULT EQUAL 0)
                                        set(HIP_GPU_AVAILABLE ON)
                                endif()
                        endif()
                endif()
                if(HIP_GPU_AVAILABLE)
                        set(BUILD_GPU_TARGET ON)
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
                        message(STATUS "No AMD GPU detected. Skipping HIP target.")
                endif()
        else()
	    message(STATUS "No HIP found")
	endif()
endif()

if (USE_CPU)
        add_executable (dipper_cpu
                        cpu/tree_generation.cpp
                        cpu/placement.cpp
                        cpu/mash.cpp
                        src/twoBitCompressor.cpp
                        src/fourBitCompressor.cpp
                        cpu/matrix_reader.cpp
                        cpu/MSA.cpp
                        src/tree.cpp
                        cpu/neighborJoining.cpp
                        cpu/placement_close_k.cpp
                        cpu/divide_and_conquer/mash_cu.cpp
                        cpu/divide_and_conquer/mash.cpp
                        cpu/divide_and_conquer/msa.cpp
                        cpu/divide_and_conquer/placement_close_k.cpp
                        cpu/divide_and_conquer/placement_close_k_cu.cpp)
        
        target_compile_options(dipper_cpu PRIVATE -w)
endif()

add_compile_definitions(TBB_PREVIEW_GLOBAL_CONTROL=1)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})
find_package(TBB REQUIRED CONFIG)

if (BUILD_GPU_TARGET)
        add_executable(dipper ${SRC_FILES})
        set_target_properties(dipper PROPERTIES CUDA_SEPARABLE_COMPILATION ON CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON) # set_target_properties(dipper_cpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
        find_library(STDCXXFS NAMES stdc++fs)
        if (STDCXXFS)
        target_link_libraries(dipper PRIVATE ${STDCXXFS}) #   target_link_libraries(dipper_cpu PRIVATE ${STDCXXFS})
        endif()
        target_link_libraries(dipper PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb) # target_link_libraries(dipper_cpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
        install(TARGETS dipper RUNTIME DESTINATION bin)
endif()

if (USE_CPU)
        set_target_properties(dipper_cpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(dipper_cpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
        install(TARGETS dipper_cpu RUNTIME DESTINATION bin)
endif()
