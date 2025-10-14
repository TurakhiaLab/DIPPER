#!/bin/bash

BUILD_TYPE=$1
CMAKE_OPTIONS=""

if [ "$BUILD_TYPE" == "cuda" ]; then
    CMAKE_OPTIONS="-DUSE_CUDA=ON -DUSE_HIP=OFF"
elif [ "$BUILD_TYPE" == "hip" ]; then
    CMAKE_OPTIONS="-DUSE_CUDA=OFF -DUSE_HIP=ON"
else
    echo "Error: either pass cuda or hip. Usage: ./installUbuntu.sh cuda/hip"
    exit 1
fi

SCRIPT_DIR=$(pwd)
BUILD_DIR="${SCRIPT_DIR}/../build"
BIN_DIR="${SCRIPT_DIR}/../"
TBB_VERSION="2021.9.0"
TBB_ARCHIVE="v${TBB_VERSION}.tar.gz"
TBB_DIR_NAME="oneTBB-${TBB_VERSION}"
TBB_INSTALL_DIR="${BUILD_DIR}/${TBB_DIR_NAME}/install"
TBB_CMAKE_DIR=""

find_tbb_cmake_dir() {
    if [[ "$(uname)" == "Darwin" ]]; then # macOS
        return 0
    fi
    for libdir in lib lib64 lib32 libx32; do
        if [ -f "${TBB_INSTALL_DIR}/${libdir}/cmake/TBB/TBBConfig.cmake" ]; then
            TBB_CMAKE_DIR="${TBB_INSTALL_DIR}/${libdir}/cmake/TBB"
            return 0
        fi
    done
    return 1
}

find_tbb_dev() {
    if dpkg -l | grep -q libtbb-dev; then
        return 0
    else
        return 1
    fi
}

install_tbb() {
    echo "Installing TBB..."
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}" || exit 1
    wget "https://github.com/oneapi-src/oneTBB/archive/refs/tags/${TBB_ARCHIVE}"
    tar -xvzf "${TBB_ARCHIVE}"
    cd "${TBB_DIR_NAME}" || exit 1
    mkdir -p build && cd build
    cmake .. -DCMAKE_CXX_FLAGS="-Wno-stringop-overflow" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="${TBB_INSTALL_DIR}"
    cmake --build . --config Release --parallel
    cmake --install .
}

if ! find_tbb_dev; then
    if ! find_tbb_cmake_dir; then
        echo "libtbb-dev not found. Install TBB from source."
        install_tbb
        if ! find_tbb_cmake_dir; then
            echo "Error: Could not find TBB CMake config after installation."
            exit 1
        fi
    else
        if [[ "$(uname)" == "Darwin" ]]; then # macOS
            if brew list --versions tbb &> /dev/null; then
                tbb_ver=$(brew list --versions tbb | cut -c5-100)
                echo "Found TBB: Version: ${tbb_ver}"
                mkdir -p "${BUILD_DIR}"
            else
                echo "TBB is not found via Homebrew. Please use "bash ./install/installDependencies.sh" to install libraries."
            fi
        else
            echo "TBB already installed at: ${TBB_CMAKE_DIR}"
        fi
    fi
else 
    tbb_ver=$(dpkg -s libtbb-dev | grep '^Version')
    echo "Found TBB: ${tbb_ver}"
    mkdir -p "${BUILD_DIR}"
fi


HIPCC_PATH=$(which hipcc)
HIP_COMPILE_VERSION=$(echo "$HIPCC_PATH" | sed -n 's|.*/rocm-\([0-9.]*\)/.*|\1|p')
echo ${HIP_COMPILE_VERSION}

cd "${BUILD_DIR}" || exit 1
rm -rf CMake*
if [[ "$(uname)" == "Darwin" ]]; then
    cmake -DTBB_DIR="$(brew --prefix tbb)/lib/cmake/tbb" -DHIP_COMPILE_VERSION=${HIP_COMPILE_VERSION} ..
else
    if find_tbb_dev; then
        cmake -DHIP_COMPILE_VERSION=${HIP_COMPILE_VERSION} $CMAKE_OPTIONS ..
    else
        cmake -DTBB_DIR="${TBB_CMAKE_DIR}" -DHIP_COMPILE_VERSION=${HIP_COMPILE_VERSION} $CMAKE_OPTIONS ..
    fi
fi
cmake --build . --config Release --parallel
# cmake --install . --prefix ${BIN_DIR}
export PATH=${BIN_DIR}/bin:$PATH
cd ${BIN_DIR}
