leaves=$1
out_dir=$2
length=$3
bl=$4
mkdir -p $out_dir


#/usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim /data/swalia/dipper/alisim/datasets/10M --length 10000 -m "GTR+G+I"  seed 1 -t "RANDOM{yh/10000000}" -nt 32 -redo -af fasta -rlen 0.0000001 0.001 0.005
# if bl is not set, use the following command 

echo "Generating alisim dataset with $leaves leaves, length $length, bl $bl"

# if [ -z "$bl" ]; then
/usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim ${out_dir}/data --length $length -m "GTR+G+I" --indel-size POW{1.7/50},POW{1.7/50}  --indel 0.03,0.09 seed 1 -t "RANDOM{yh/${leaves}}" -nt 1 -rlen 0.000002 0.00002 0.0002;
   # -t /data/swalia/dipper/alisim/trees/paper/nleaves${leaves}.treefile -nt 32;
# else
#     /usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim ${out_dir}/len_${length}_leaves_${leaves} --length $length -m "GTR+G+I" --indel-size POW{1.7/50},POW{1.7/50}  --indel 0.03,0.09 seed 1 -t /data/swalia/dipper/alisim/trees/paper/bl_effect/nleaves${leaves}_bl_${bl}.treefile -nt 32;
# fi
#/usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim ${out_dir}/len_${length}_leaves_${leaves} --length $length -m "GTR" --indel-size POW{1.7/50},POW{1.7/50} --indel 0.03,0.09 seed 1 -t /data/swalia/dipper/ALISIM/benchmarks/input/normalbranches_nLeaves${leaves}.treefile -nt 32;
#/usr/bin/time -v timeout 1d /home/swalia@AD.UCSD.EDU/tools/iqtree-2.3.5-Linux-intel/bin/iqtree2 --alisim ${out_dir}/len_${length}_leaves_${leaves} --length $length -m "GTR{0.013206908228919744,0.10497798848628515,0.04165255672197765,0.007450050795800881,1.0253979004402303}+F{0.3,0.2,0.2,0.3}+G4{0.5}+I{0.2}" --indel 0.03,0.09 seed 1 -t /data/swalia/dipper/ALISIM/benchmarks/input/normalbranches_nLeaves${leaves}.treefile -nt 32;
