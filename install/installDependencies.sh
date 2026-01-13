#!/bin/bash

OS_NAME="$(uname -s)"

if [[ "${OS_NAME}" == "Darwin" ]]; then
    if ! command -v brew >/dev/null 2>&1; then
        echo "Homebrew not found. Install Homebrew from https://brew.sh/"
        exit 1
    fi
    brew update
    brew install tbb boost cmake wget git
else
    sudo apt-get update && \
        sudo apt-get install -y \
            wget \
            git \
            vim \
            build-essential \
            libboost-all-dev \
            libtbb-dev \
            cmake && \
        sudo apt-get clean
fi
