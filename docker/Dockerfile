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
ENV CONDA_DIR /opt/conda

# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH

# Set up Conda channels and create environment
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r    && \
    conda create -n dipper python=3.11 -y
RUN echo "source activate dipper" > ~/.bashrc
SHELL ["conda", "run", "-n", "dipper", "/bin/bash", "-c"]

# Clone and install DIPPER inside the dipper environment
RUN git clone https://github.com/TurakhiaLab/DIPPER.git 
RUN conda run -n dipper /bin/bash -c "cd /home/DIPPER/install && ./installUbuntu.sh cuda"

ENV PATH="/opt/DIPPER/bin:$PATH"

RUN conda run -n dipper conda install -y rapidnj fastme decenttree quicktree ccphylo fastphylo && \
    conda run -n dipper pip install apples

WORKDIR /home/DIPPER/bin

CMD ["/bin/bash"]