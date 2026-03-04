#!/bin/bash

PROJECT="dom_admixed_rm"
FILE="./out/s02.climate_data/orchard/climate_extract_current_orchard.csv"
REPEAT_NUMBER=89
OUT_FOLD="./out/s03.orchard_generalization/${PROJECT}"

mkdir -p "${OUT_FOLD}"

# Using sed to remove quotes and skip the header, then processing with while
# We change IFS to a comma to handle the CSV format
sed '1d; s/"//g' "$FILE" | while IFS=',' read -r Site bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19; do

    # Define the output file path
    output_file="${OUT_FOLD}/env_${Site}_Chelsa"
    
    # Create/Clear the file
    : > "$output_file"

    # Prepare the data line (tab-separated as per your echo -e)
    data_line=$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
        "$bio1" "$bio2" "$bio3" "$bio4" "$bio5" "$bio6" "$bio7" "$bio8" "$bio9" \
        "$bio10" "$bio11" "$bio12" "$bio13" "$bio14" "$bio15" "$bio16" "$bio17" "$bio18" "$bio19")

    # Efficiently append the lines
    for i in $(seq 1 $REPEAT_NUMBER); do
        echo -e "$data_line" >> "$output_file"
    done

    echo "Generated: env_${Site}_Chelsa"
done