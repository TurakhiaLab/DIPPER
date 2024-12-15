# phylo-accel

## Make sure going to placement subfolder. Everything is there.

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make mash-placement
```

## Run Instructions
```
./mash-placement -f [input-file] -i [r or m or d] -o t > [output-file]
# -i r -> raw sequences
# -i m -> MSA
# -i d -> distance matrix
# for more options, try ./mash-placement -h
```

## Outputs

```
Number -x occurs at indices: a b h t 
Number -y occurs at indices: o i j p
```


where x and y are cluster numbers, a b h t are sequences belonging to cluster x, and similarly o i j p belong to cluster y. 

The output for each cluster's tree is also printed in the Newick format. 



The build instruction is the same as the original. So the following statements are going to be unchanged. After running the workflow multiple sub-trees/clusters will be formed and printed out in the following manner 

## Example Input and Output
### Input
```
>Seq1
AGTACAATGACATCTCGCTCTGCCCCCGTCCAGACTGCCCACAATTGGAC
>Seq2
AGTACAATGACATCTCGCTCTGCCCTCATCCAGAGTGGCCACAATTCGAC
>Seq3
GGTACAATAACATCTCTCTCTGCCCCTGTCCAAACTGCCTACAATTGGAC
>Seq4
AGTACAATGACATCACGATCTGCCCCCGTCCAGGATGCCCACCATTGGAC
>Seq5
AGTACATTGAAATCTCGCTCTGCCCCCGTCCAGACCGCCCACAATTGGAC
>Seq6
AGTACAATGACATCTCACTCTGCCCCCGTCCAGACTGCCCACAATTGGAC
>Seq7
AGTACAATGCCATCTCGCTCTGCCCCCGTCCAGACTGCCCACAATTGGAC
>Seq8
AGTACAATGACATCTCGCTCTCCCCCCCTCCAGACTGCCCACAATTGGAC
>Seq9
AGTACAATGAAATCTCGATCTGCCCCCGTCCAGACTGCCCACAATTGGAC
>Seq10
AGTACAATGACATCTCGCTCAGCCCCCGTCCAGAATTCCCGAAATTGGAC

```

### Output
```
Unique numbers and their indices:
Number -2 occurs at indices: 4 5 7 9 
Number -1 occurs at indices: 0 1 2 3 6 8 

Seq8 Seq2 Seq1 Seq7 Distance Operation Time 0 ms
Tree Operation Time 0 ms
((Seq7:0.00061659,Seq2:0.00058156):1.7511e-05,(Seq1:0.00025149,Seq8:0.00035672):0.00024235);
Tree Created in: 1 ms

Seq6 Seq4 Seq9 Seq5 Seq10 Seq3 Distance Operation Time 0 ms
Tree Operation Time 1 ms
((Seq5:0.00061659,Seq4:0.00058156):1.7511e-05,((Seq3:0.00063407,Seq10:0.00035411):0.00019279,(Seq9:0.00025149,Seq6:0.00035672):0.00015499):8.7368e-05);
Tree Created in: 1 ms
```

