#!/bin/bash

# Input file
input_file="./out/Bioclim_orchard.csv"

# Number of times to repeat columns 2 to 20
repeat_count=64  # Adjust this to the number of repetitions you want

# Initialize a counter to skip the first line
header_skipped=false

# Read the input file line by line
while IFS=$'\t' read -r site bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19; do
    # Skip the header line
    if [ "$header_skipped" = false ]; then
        header_skipped=true
        continue
    fi

    # Create a file named after the site
    output_file="./out/env_${site}_Chelsa"
    > "$output_file"  # Clear or create the output file

    # Append repeated columns with a newline after each repetition
    for i in $(seq 1 $repeat_count); do
        echo -e "$bio1\t$bio2\t$bio3\t$bio4\t$bio5\t$bio6\t$bio7\t$bio8\t$bio9\t$bio10\t$bio11\t$bio12\t$bio13\t$bio14\t$bio15\t$bio16\t$bio17\t$bio18\t$bio19" >> "$output_file"
    done

done < "$input_file"

