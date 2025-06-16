#!/bin/bash

SCRIPT_DIR=$(pwd)
BUILD_DIR="${SCRIPT_DIR}/../build"
BIN_DIR="${SCRIPT_DIR}/../bin"

mkdir -p ${BUILD_DIR}
cd "${BUILD_DIR}" || exit 1
cmake ..
cmake --build . --config Release --parallel
cmake --install . --prefix ${BIN_DIR}
export PATH=${BIN_DIR}:$PATH