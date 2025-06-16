FROM nvidia/cuda:12.2.2-devel-ubuntu22.04

ENV DEBIAN_FRONTEND noninteractive
RUN chmod 1777 /tmp

WORKDIR /home/

# install dependencies
RUN apt-get update && \
    apt-get install -y \
       	wget \
	    git \
        vim \
	    build-essential \
	    libboost-all-dev \
		libtbb-dev \
	    cmake && \
	apt-get clean

# Install DIPPER
RUN git clone https://github.com/TurakhiaLab/DIPPER.git && \
    cd DIPPER/scripts && \
    chmod +x installUbuntu.sh && \
    ./installUbuntu.sh && \
    export PATH=/home/DIPPER/bin/:$PATH





