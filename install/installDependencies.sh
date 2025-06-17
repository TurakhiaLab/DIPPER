#!/bin/bash

sudo    apt-get update && \
        apt-get install -y \
            wget \
            git \
            vim \
            build-essential \
            libboost-all-dev \
            libtbb-dev \
            cmake && \
        apt-get clean