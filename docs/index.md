# <b>Welcome to DIPPER Wiki</b>
<div align="center">
    <img src="images/logo.png"/ width="300">
</div>

## <b>Introduction</b> 
### <b>Overview</b><a name="overview"></a>

DIPPER (**DI**stance-based **P**hylogenetic **P**lac**ER**) is a tool for ultrafast and ultralarge phylogenetic reconstruction on GPUs, designed to maintain high accuracy with a minimal memory footprint. DIPPER introduces several innovations, including a divide-and-conquer strategy, a new placement algorithm, and an on-the-fly distance calculator that dynamically enables selective distance computation. DIPPER consistently outperforms existing distance-based methods in speed, accuracy, and memory efficiency. In addition, DIPPER minimizes branch length underestimation for non-additive distance matrices compared to earlier methods and offers a strict mode that completely eliminates the underestimation. 

<div align="center">
    <div><b>Figure 1: Overview of DIPPER algorithm</b></div>
    <img src="images/overview.png" width="1000"/>
</div>


<a name="install"></a>
## <b>Installation Methods</b>

NOTE: DIPPER is currently supported on systems with <b>NVIDIA GPUs only</b>. Support for additional platforms, including AMD GPUs and CPU-only options for x86-64 and ARM64 architecture, will be added soon. Stay tuned!

### 1. <a name="conda"></a> Using Conda
DIPPER is available on platforms with NVIDIA GPUs via Conda. See [DIPPER Bioconda Page](https://anaconda.org/bioconda/dipper) for details.
#### i. Dependencies
1. [Conda](https://docs.conda.io/en/latest/)

#### ii. Create and activate a Conda environment
```bash
conda create -n dipper python=3.11 -y
conda activate dipper
# Set up channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
# Install DIPPER
conda install bioconda::dipper
```

#### iii. Run DIPPER
```bash
# Insider conda environment
dipper --help
```

### 2. <a name="dockerimage"></a> Using Docker Image
To use DIPPER in a docker container, users can create a docker container from a docker image, by following these steps
#### i. Dependencies
1. [Docker](https://docs.docker.com/engine/install/)
#### ii. Pull and build the DIPPER docker image from DockerHub
```bash
## Note: If the Docker image already exists locally, make sure to pull the latest version using 
## docker pull swalia14/dipper:latest

## If the Docker image does not exist locally, the following command will pull and run the latest version
docker run -it --gpus all swalia14/dipper:latest
```
#### iii. Run DIPPER
```bash
# Insider docker container (path: /home/Dipper/bin)
./dipper --help
```

### 3. Using DockerFile <a name="dockerfile"></a>
Docker container with the preinstalled DIPPER program can also be built from a Dockerfile by following these steps.

#### i. Dependencies
1. [Docker](https://docs.docker.com/engine/install/)
2. [Git](https://git-scm.com/downloads)

#### ii. Clone the repository and build a docker image
```bash
git clone https://github.com/TurakhiaLab/DIPPER.git
cd DIPPER/docker
docker build -t dipper .
```
#### iii. Build and run the docker container
```bash
docker run -it --gpus all dipper
```
#### iv. Run DIPPER
```bash
# Insider docker container (path: /home/Dipper/bin)
./dipper --help
```

### 4. <a name="script"></a> Using installation script (requires sudo access)  

Users without sudo access are advised to install DIPPER via [Docker Image](#dockerimage) or [Dockerfile](#dockerfile).

**Step 1:** Clone the repository
```bash
git clone https://github.com/TurakhiaLab/DIPPER.git
cd DIPPER
```
**Step 2:** Install dependencies (requires sudo access)

DIPPER depends on the following common system libraries, which are typically pre-installed on most development environments:
```bash
- wget
- cmake 
- build-essential 
- libboost-all-dev
- libtbb-dev
```
For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt install -y wget cmake build-essential libboost-all-dev  libtbb-dev
```

**Step 3:** Build DIPPER

```bash
cd install
chmod +x installUbuntu.sh
./installUbuntu.sh
cd ../
```

**Step 4:** The DIPPER executable is located in the `bin` directory and can be run as follows:
```bash
cd bin
./dipper --help
```

## <a name="Run DIPPER"></a> **Run DIPPER**

### Functionalities

<!-- <div align="center"> -->
<div name="table1"> <b>Table 1:</b> List of functionalities supported by DIPPER </div>

| **Option**               | **Description**                                                                                                  |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------- |
| `-i`, `--input-format`   | Input format <b>(required)</b>:<br>`d` - distance matrix<br>`r` - raw sequences<br>`m` - MSA                                       |
| `-o`, `--output-format`  | Output format:<br>`t` - phylogenetic tree in Newick format (default) <br>`d` - distance matrix in PHYLIP format (coming soon!!)                                              |
| `-I`, `--input-file`     | Input file path <b>(required)</b>:<br>PHYLIP for distance matrix, FASTA for MSA or raw sequences                                   |
| `-O`, `--output-file`    | Output file path <b>(required)</b>                                                                                                |
| `-m`, `--algorithm`      | Algorithm selection:<br>`0` - auto <b>(default)</b><br>`1` - force placement<br>`2` - force NJ<br>`3` - divide-and-conquer     |
| `-K`, `--K-closest` | Placement mode:<br>`-1` - exact<br>`10` - <b>default</b>                                                      |
| `-k`, `--kmer-size`      | K-mer size (Valid range: 2–15, <b>default: 15</b>)                                                                      |
| `-s`, `--sketch-size`    | Sketch size (<b>default: 1000</b>)                                                                                    |
| `-d`, `--distance-type`  | Distance type:<br>`1` - uncorrected<br>`2` - JC <b>(default)</b> <br>`3` - Tajima-Nei<br>`4` - K2P<br>`5` - Tamura<br>`6` - Jinnei |
| `-a`, `--add`            | Add query to a backbone tree using k-closest placement                                                           |
| `-t`, `--input-tree`     | Input backbone tree in Newick <b>(required with --add option)</b> format                                                                             |
| `-h`, `--help`           | Show help message                                                                                                |

<!-- </div> -->

!!!Note
    All the files in the examples below can be found in the `DIPPER/dataset`.

Enter into the build directory (assuming `$DIPPER_HOME` directs to the DIPPER repository directory)  
```bash
cd $DIPPER_HOME/bin
./dipper -h
```

### De-novo phylogeny construction <a name="denovo"></a>
DIPPER supports de-novo construction of phylogenies from unaligned/aligned sequences in FASTA format and distance matrix in PHYLIP format. 

#### Default mode <a name="default"></a>
In default mode, DIPPER constructs phylogeny using:
1. Conventional NJ for sequences/tips < 30,000
2. Placement technique for sequences/tips >= 30,000 and < 1,000,000
3. Divide-and-conquer technique for sequences/tips >= 1,000,000

##### From unaligned sequences
Usage syntax
```bash
./dipper -i r -o t -I <path to unaligned sequences FASTA file> -O <path to output file>
```
Example
```bash
./dipper -i r -o t -I ../dataset/t2.unaligned.fa -O tree.nwk
```

##### From aligned sequences
Usage syntax (using JC model)
```bash
./dipper -i m -o t -d 2 -I <path to aligned sequences FASTA file> -O <path to output file>
```
Example
```bash
./dipper -i m -o t -d 2 -I ../dataset/t1.aligned.fa -O tree.nwk
```

##### From distance matrix
Usage syntax 
```bash
./dipper -i d -o t -I <path to distance matrix PHYLIP file> -O <path to output file>
```
Example
```bash
./dipper -i d -o t -I ../dataset/t2.phy -O tree.nwk
```

#### Construct phylogeny using placement technique <a name="place"></a>
DIPPER allows users to construct phylogeny using the forced placement technique by setting the `-m` option to `1`. Below we provide a syntax and an example for input unaligned sequences, but DIPPER also supports aligned sequences and distance matrix as input.
Usage syntax
```bash
./dipper -i r -o t -m 1 -I <path to unaligned sequences FASTA file> -O <path to output file>
```
Example
```bash
./dipper -i r -o t -m 1 -I ../dataset/t2.unaligned.fa -O tree.nwk
```

#### Construct phylogeny using divide-and-conquer technique <a name="dc"></a>
DIPPER allows users to construct phylogeny using the forced divide-and-conquer technique by setting the `-m` option to `3`. Below we provide a syntax and an example for input unaligned sequences, but DIPPER also supports aligned sequences and distance matrix as input.
Usage syntax
```bash
./dipper -i r -o t -m 3 -I <path to unaligned sequences FASTA file> -O <path to output file>
```
Example
```bash
./dipper -i r -o t -m 3 -I ../dataset/t2.unaligned.fa -O tree.nwk
```

### Adding tips (sequences) to a backbone tree <a name="add"></a>
DIPPER allows users to add tips to an existing backbone tree using the placement technique. It requires tip sequences from the backbone tree and input query sequences to be provided in a single file (FASTA format), along with the input tree in Newick format.

Usage syntax 
```bash
./dipper -i r -o t -m 1 --add -I <path to unaligned/aligned sequences FASTA file (containing backbone tree tip sequences and query sequences)> -O <path to output file> -t <path to input tree>
```
Example
```bash
./dipper -i r -o t -m 1 --add -I ../dataset/t2.unaligned.fa -O tree.nwk -t ../dataset/backbone.nwk
```


##  <a name="contribution"></a> Contributions
We welcome contributions from the community to enhance the capabilities of **DIPPER**. If you encounter any issues or have suggestions for improvement, please open an issue on [DIPPER GitHub page](https://github.com/TurakhiaLab/DIPPER/issues). For general inquiries and support, reach out to our team.

##  <a name="cite"></a> Citing DIPPER
If you use DIPPER in your research or publications, we kindly request that you cite the following paper:         
* Sumit Walia, Zexing Chen, Yu-Hsiang Tseng, Yatish Turakhia, "<i>Ultrafast and Ultralarge Distance-Based Phylogenetics Using DIPPER</i>", bioRxiv 2025.08.12.669583; doi: [https://doi.org/10.1101/2025.08.12.669583](https://doi.org/10.1101/2025.08.12.669583)
