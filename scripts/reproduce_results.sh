# This script reproduces the results described here: https://zenodo.org/records/17259722

# Download datasets
wget https://zenodo.org/records/17259722/files/DIPPER_data.tar.gz?download=1 -O DIPPER_data.tar.gz
tar -xvzf DIPPER_data.tar.gz
rm DIPPER_data.tar.gz

############################### AliSim dataset (unaligned sequences) ######################################

## Directory structure 
## AliSim
## └── dataset
##     ├── 10K
##     │    ├── data.unaligned.fasta.gz
##     │    └── data.treefile.nwk.gz (reference tree)
##     ├── 20K
##     ├── ..
##     ├── 1M
##     └── ultralarge
##         ├── 200K
##         ├── 500K
##         ├── ..
##         └── 10M

./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/10K/data.unaligned.fa.gz -O dipper_alisim_10K -m 2
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/20K/data.unaligned.fa.gz -O dipper_alisim_20K -m 2
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/50K/data.unaligned.fa.gz -O dipper_alisim_50K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/100K/data.unaligned.fa.gz -O dipper_alisim_100K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/200K/data.unaligned.fa.gz -O dipper_alisim_200K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/500K/data.unaligned.fa.gz -O dipper_alisim_500K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/1M/data.unaligned.fa.gz -O dipper_alisim_1M -m 1

### The ultralarge datasets
#### Placement mode
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/200K/data.unaligned.fa.gz -O dipper_alisim_ultralarge_200K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/500K/data.unaligned.fa.gz -O dipper_alisim_ultralarge_500K -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/1M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_1M -m 1
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/2M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_2M -m 1
#### Divide-and-conquer mode
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/200K/data.unaligned.fa.gz -O dipper_alisim_ultralarge_200K_dac -m 3
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/500K/data.unaligned.fa.gz -O dipper_alisim_ultralarge_500K_dac -m 3
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/1M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_1M_dac -m 3
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/2M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_2M_dac -m 3
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/5M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_5M_dac -m 3
./dipper -i r -I $PWD/DIPPER_data/AliSim/dataset/ultralarge/10M/data.unaligned.fa.gz -O dipper_alisim_ultralarge_10M_dac -m 3


############################### RNASim dataset (aligned sequences) ######################################
    
## Directory structure
## RNASim
## └── dataset
##     ├── 10K
##     │    ├── data.aligned.fasta.gz
##     │    └── data.treefile.nwk.gz (reference tree)
##     ├── 20K
##     ├── ..
##     └── 1M

./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/10K/data.aligned.fa.gz -O dipper_rnasim_10K -m 2
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/20K/data.aligned.fa.gz -O dipper_rnasim_20K -m 2
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/50K/data.aligned.fa.gz -O dipper_rnasim_50K -m 1
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/100K/data.aligned.fa.gz -O dipper_rnasim_100K -m 1
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/200K/data.aligned.fa.gz -O dipper_rnasim_200K -m 1
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/500K/data.aligned.fa.gz -O dipper_rnasim_500K -m 1
./dipper -i m -I $PWD/DIPPER_data/RNASim/dataset/1M/data.aligned.fa.gz -O dipper_rnasim_1M -m 1


################################## SILVA datasets ##################################################

### lsu subunit
### Download sequences
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_tax_silva_full_align_trunc.fasta.gz -O SILVA_138.2_SSURef_tax_silva_full_align_trunc.fasta.gz
./dipper -i m -I $PWD/DIPPER_data/SILVA/SILVA_138.2_SSURef_tax_silva_full_align_trunc.fasta.gz -O dipper_silva_lsu -m 1

### ssu subunit
### Download sequences
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva_full_align_trunc.fasta.gz -O SILVA_138.2_LSURef_NR99_tax_silva_full_align_trunc.fasta.gz
./dipper -i m -I $PWD/DIPPER_data/SILVA/SILVA_138.2_LSURef_NR99_tax_silva_full_align_trunc.fasta.gz -O dipper_silva_ssu -m 3