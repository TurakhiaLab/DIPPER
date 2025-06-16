#!/bin/bash

startDir=$(pwd)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BUILD_DIR="${SCRIPT_DIR}/../build"

mkdir -p ${BUILD_DIR}
cd "${BUILD_DIR}" || exit 1
cmake ..
