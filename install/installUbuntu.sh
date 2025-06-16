#!/bin/bash

SCRIPT_DIR=$(pwd)
BUILD_DIR="${SCRIPT_DIR}/../build"

mkdir -p ${BUILD_DIR}
cd "${BUILD_DIR}" || exit 1
cmake ..
