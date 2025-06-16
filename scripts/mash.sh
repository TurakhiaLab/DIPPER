#!/bin/bash
# Ex: source mash.sh mash len_1000.dist
mash_directory=$1
out_file=$2
files=("$mash_directory"/*.fa)
num_files=${#files[@]}
declare -A dist_matrix
kmer=15
THREADS=32

for ((i=0; i<num_files; i++)); do
    a=${files[i]}
    base_a=$(basename "$a" .fa)
    mash sketch -p $THREADS "$a" -o "$base_a" -k $kmer
done

for ((i=0; i<num_files; i++)); do
    a=${files[i]}
    base_a=$(basename "$a" .fa)
    printf "%s\n" "$i"
    dist_matrix[$i,$i]=0.00
    for ((j=i+1; j<num_files; j++)); do
        b=${files[j]}
        base_b=$(basename "$b" .fa)
        dist=$(mash dist -p $THREADS "$base_a.msh" "$base_b.msh" | cut -f3)    
        dist_matrix[$i,$j]=$dist
    done
done



# for ((i=0; i<num_files; i+=$THREADS)); do
#     for ((k=i; k<i+$THREADS && k<num_files; k++)); do
#         a=${files[k]}
#         base_a=$(basename "$a" .fa)
#         {
#             mash sketch "$a" -o "$base_a" -k $kmer
#         } &
#     done
#     wait
# done

# for ((i=0; i<1; i++)); do
#     a=${files[i]}
#     base_a=$(basename "$a" .fa)
#     printf "%s\n" "$i"
#     dist_matrix[$i,$i]=0.00
#     for ((j=i+1; j<3; j+=2)); do
#         dist_matrix[$i,$j]=$(mash dist "$base_a.msh" "$(basename "${files[j]}" .fa).msh" | cut -f3) &
#         dist_matrix[$i,$((j+1))]=$(mash dist "$base_a.msh" "$(basename "${files[j+1]}" .fa).msh" | cut -f3) &
#         wait
#         echo ${dist_matrix[$i,$j]}
#         echo ${dist_matrix[$i,$((j+1))]}
#     done
# done

rm *.msh

# Print the 2D array
string="$num_files\n"
for ((i=0; i<num_files; i++)); do
    a=${files[i]}
    base_a=$(basename "$a" .fa)
    string+="$base_a\t"
    for ((j=0; j<num_files; j++)); do
        if [[ $j -lt $i ]]; then
            string+="${dist_matrix[$j,$i]}\t"
        else 
            string+="${dist_matrix[$i,$j]}\t"
        fi
    done
    string+="\n"
    echo -e "$string" >> $out_file
    string=""
done
