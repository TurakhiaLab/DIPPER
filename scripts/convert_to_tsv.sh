
printf "\t"
cut -f1 $1 | tr "\n" "\t" | sed "s/\t$/\n/g"

#printf "%s\t" "${1%.*}" | sed "s/.*dist\-//g" 
#printf "%s\t" "mix" 
cut -f2 $1 | tr "\n" "\t" | sed "s/\t$/\n/g"

