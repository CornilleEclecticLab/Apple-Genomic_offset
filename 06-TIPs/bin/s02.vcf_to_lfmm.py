import sys
import gzip
import shutil

def convert_genotype(gt):
    if gt in {"./.", ".|."}:
        return "9"   # PLINK missing
    sep = "/" if "/" in gt else "|"
    alleles = gt.split(sep)
    try:
        alleles = list(map(int, alleles))
    except ValueError:
        return "9"
    return str(sum(alleles))

def determine_variant_type(var_id):
    if ".INS." in var_id:
        return "INS"
    elif ".DEL." in var_id:
        return "DEL"
    else:
        return "OTHER"

def parse_vcf(input_path, prefix):
    open_func = gzip.open if input_path.endswith(".gz") else open

    variants = []
    genotypes_matrix = []
    samples = []

    with open_func(input_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                samples = header[9:]
                continue

            fields = line.strip().split("\t")
            chrom, pos, var_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            variant_type = determine_variant_type(var_id)
            if variant_type not in {"INS", "DEL"}:
                print(var_id)
                continue  # keep only indels

            # Save variant info
            variants.append((chrom, var_id))

            # Genotypes for all samples
            genotypes = [field.split(":")[0] for field in fields[9:]]
            converted = [convert_genotype(gt) for gt in genotypes]
            genotypes_matrix.append(converted)

    # Transpose matrix: variants x samples → samples x variants
    transposed = list(zip(*genotypes_matrix))

    # Write .012
    with open(prefix + ".012", "w") as fout:
        for row in transposed:
            fout.write(" ".join(row) + "\n")

    # Write .indv
    with open(prefix + ".012.indv", "w") as fout:
        for s in samples:
            fout.write(s + "\n")

    # Write .pos
    with open(prefix + ".012.pos", "w") as fout:
        for chrom, var_id in variants:
            fout.write(f"{chrom}\t{var_id}\n")

    # Make copy as .lfmm
    out_012 = prefix + ".012"
    out_lfmm = prefix + ".lfmm"
    shutil.copy(out_012, out_lfmm)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vcf_to_lfmm.py input.vcf[.gz] output_prefix")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_prefix = sys.argv[2]
    parse_vcf(input_vcf, output_prefix)
