#!/usr/bin/env bash

module load all gencore/3
module load vcftools

##############################################################
# USER CONFIGURATION
##############################################################

Project="tips_dom_admixed_rm_V2" # Your workspace
VCF="./out/${Project}/s01.filter_vcf/GraffiTE.merged.genotypes_${Project}.vcf.gz" # Your TIPs VCF file
POP_TABLE="./in/group_geo_tips_dom_admixed_rm.txt"  # The file you pasted
POP_DIR="./out/$Project/popgroups"             # Output folder for group files

mapfile -t GROUPS < <(for f in "$POP_DIR"/*; do basename "$f"; done)

# show them and count
echo "Detected populations:"
printf '%s\n' "${GROUPS[@]}"
echo "Number of populations: ${#GROUPS[@]}"

# ##############################################################
# # PREP
# ##############################################################

# echo "== Creating population directory =="
# mkdir -p $POP_DIR

# echo "== Parsing population table and generating group files =="

# # Reset population directory
# rm -f $POP_DIR/*

# # Generate one file per population
# awk '
# {
#     ind=$1;
#     pop=$2;
#     out=sprintf("%s/%s.txt","'"$POP_DIR"'",pop);
#     print ind >> out;
# }
# ' $POP_TABLE

# echo "== Done generating group files =="
# echo

# # List of groups
# #GROUPS=($(ls $POP_DIR))
# mapfile -t GROUPS < <(ls "$POP_DIR")

# echo "Detected populations:"
# printf '%s\n' "${GROUPS[@]}"
# echo

##############################################################
# OUTPUT MATRIX FILES
##############################################################

FST_MATRIX="./out/$Project/s08.tips_fst_dxy/FST_matrix.tsv"
DXY_MATRIX="./out/$Project/s08.tips_fst_dxy/dXY_matrix.tsv"

# echo -e "group\t${GROUPS[*]}" > $FST_MATRIX
# echo -e "group\t${GROUPS[*]}" > $DXY_MATRIX

# ##############################################################
# # FUNCTIONS
# ##############################################################

# compute_dxy () {
#     local pop1=$1
#     local pop2=$2

#     # Compute allele freq for each group
#     vcftools --vcf $VCF --keep $POP_DIR/$pop1 --freq --out freq_$pop1 >/dev/null
#     vcftools --vcf $VCF --keep $POP_DIR/$pop2 --freq --out freq_$pop2 >/dev/null

#     # Join and calculate mean dXY
#     DXY=$(join \
#         <(awk 'NR>1 && $5!~"nan"{print $1":"$2, $5}' freq_${pop1}.frq | sort) \
#         <(awk 'NR>1 && $5!~"nan"{print $1":"$2, $5}' freq_${pop2}.frq | sort) \
#         | awk '{p1=$2; p2=$3; d=p1*(1-p2)+p2*(1-p1); s+=d; n++} END{if(n>0) print s/n; else print "NA"}' )

#     echo "$DXY"
# }

# compute_fst () {
#     local pop1=$1
#     local pop2=$2

#     vcftools \
#       --vcf $VCF \
#       --weir-fst-pop $POP_DIR/$pop1 \
#       --weir-fst-pop $POP_DIR/$pop2 \
#       --out fst_tmp >/dev/null

#     # Mean FST
#     FST=$(awk 'NR>1 && $3!="nan"{s+=$3; n++} END{if(n>0) print s/n; else print "NA"}' fst_tmp.weir.fst)

#     echo "$FST"
# }

# ##############################################################
# # MAIN LOOP
# ##############################################################

# for g1 in ${GROUPS[@]}; do

#     FST_LINE="$g1"
#     DXY_LINE="$g1"

#     for g2 in ${GROUPS[@]}; do
#         echo "Comparing $g1 vs $g2"

#         if [ "$g1" == "$g2" ]; then
#             FST_LINE="$FST_LINE\t0"
#             DXY_LINE="$DXY_LINE\t0"
#         else
#             FST_VAL=$(compute_fst $g1 $g2)
#             DXY_VAL=$(compute_dxy $g1 $g2)
#             FST_LINE="$FST_LINE\t$FST_VAL"
#             DXY_LINE="$DXY_LINE\t$DXY_VAL"
#         fi
#     done

#     echo -e "$FST_LINE" >> $FST_MATRIX
#     echo -e "$DXY_LINE" >> $DXY_MATRIX
# done

# ##############################################################
# # CLEANUP
# ##############################################################

# rm -f freq_*.frq fst_tmp.*

# echo
# echo "=================================================="
# echo "Finished!"
# echo "FST matrix: $FST_MATRIX"
# echo "dXY matrix: $DXY_MATRIX"
# echo "=================================================="


# ---------------------------
# Basic checks
# ---------------------------
if ! command -v vcftools >/dev/null 2>&1; then
  echo "ERROR: vcftools not found in PATH. Install vcftools first." >&2
  exit 1
fi

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF file '$VCF' not found." >&2
  exit 1
fi

if [ ! -f "$POP_TABLE" ]; then
  echo "ERROR: population table '$POP_TABLE' not found." >&2
  exit 1
fi

# ---------------------------
# Create pop files from table
# ---------------------------
mkdir -p "$POP_DIR"
# remove old files safely
find "$POP_DIR" -maxdepth 1 -type f -print0 | xargs -0 -r rm --

awk -v outdir="$POP_DIR" 'BEGIN{FS="[ \t]+"}
NF>=2{
  ind=$1; pop=$2;
  # print sample id (first field) into file named after population
  print ind >> (outdir"/"pop".txt")
}
' "$POP_TABLE"

# ---------------------------
# Read population filenames safely
# ---------------------------
# Use find to list files and mapfile to store them as literal strings
mapfile -t POP_FILES < <(find "$POP_DIR" -maxdepth 1 -type f -printf '%f\n' | sort)

# Remove extension (.txt) to get group names (handles files with dots robustly)
GROUPS=()
for f in "${POP_FILES[@]}"; do
  name="${f%.*}"
  GROUPS+=("$name")
done

if [ ${#GROUPS[@]} -eq 0 ]; then
  echo "ERROR: no population files found in $POP_DIR" >&2
  exit 1
fi

echo "Detected populations (${#GROUPS[@]}):"
printf '  %s\n' "${GROUPS[@]}"
echo

# ---------------------------
# Initialize output matrices
# ---------------------------
# Build header: group <tab> g1 <tab> g2 ...
header="group"
for g in "${GROUPS[@]}"; do
  header+=$'\t'"$g"
done
printf '%s\n' "$header" > "$FST_MATRIX"
printf '%s\n' "$header" > "$DXY_MATRIX"

# ---------------------------
# Helper: create temp dir
# ---------------------------
TMPDIR="$(mktemp -d)"
cleanup() {
  rm -rf "$TMPDIR"
}
trap cleanup EXIT

# ---------------------------
# Functions
# ---------------------------
compute_dxy() {
  local pop1="$1"
  local pop2="$2"
  local out1="${TMPDIR}/freq_${pop1}.frq"
  local out2="${TMPDIR}/freq_${pop2}.frq"

  vcftools --vcf "$VCF" --keep "$POP_DIR/${pop1}.txt" --freq --out "${TMPDIR}/freq_${pop1}" >/dev/null 2>&1
  vcftools --vcf "$VCF" --keep "$POP_DIR/${pop2}.txt" --freq --out "${TMPDIR}/freq_${pop2}" >/dev/null 2>&1

  # compute mean dXY; only use sites present in both freq outputs and skip nan
  dxy_val=$(join \
    <(awk 'NR>1 && $5!~"nan"{print $1":"$2, $5}' "$out1" | sort) \
    <(awk 'NR>1 && $5!~"nan"{print $1":"$2, $5}' "$out2" | sort) \
    | awk '{p1=$2; p2=$3; d=p1*(1-p2)+p2*(1-p1); s+=d; n++} END{if(n>0) printf("%.6f", s/n); else print "NA"}' )
  printf '%s' "$dxy_val"
}

compute_fst() {
  local pop1="$1"
  local pop2="$2"
  local fst_out="${TMPDIR}/fst_tmp.weir.fst"

  vcftools --vcf "$VCF" \
    --weir-fst-pop "$POP_DIR/${pop1}.txt" \
    --weir-fst-pop "$POP_DIR/${pop2}.txt" \
    --out "${TMPDIR}/fst_tmp" >/dev/null 2>&1

  # compute mean weighted FST across sites, skip nan
  fst_val=$(awk 'NR>1 && $3!="nan"{s+=$3; n++} END{if(n>0) printf("%.6f", s/n); else print "NA"}' "$fst_out")
  printf '%s' "$fst_val"
}

# ---------------------------
# Main pairwise loop
# ---------------------------
for g1 in "${GROUPS[@]}"; do
  FST_LINE="$g1"
  DXY_LINE="$g1"

  for g2 in "${GROUPS[@]}"; do
    echo "Processing: $g1 vs $g2" >&2

    if [ "$g1" = "$g2" ]; then
      FST_LINE+=$'\t'"0"
      DXY_LINE+=$'\t'"0"
    else
      fst_v=$(compute_fst "$g1" "$g2")
      dxy_v=$(compute_dxy "$g1" "$g2")
      FST_LINE+=$'\t'"$fst_v"
      DXY_LINE+=$'\t'"$dxy_v"
    fi
  done

  printf '%s\n' "$FST_LINE" >> "$FST_MATRIX"
  printf '%s\n' "$DXY_LINE" >> "$DXY_MATRIX"
done

echo
echo "Done. Outputs:"
echo "  - $FST_MATRIX"
echo "  - $DXY_MATRIX"