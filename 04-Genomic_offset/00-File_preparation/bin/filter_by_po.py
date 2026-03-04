import gzip
import sys

def LD_prunning_filtering(input_vcf_gz, prune_in_file, output_vcf):
    # Read the prune.in file into a set
    pruned_snps = set()
    with open(prune_in_file, 'r') as prune_file:
        for line in prune_file:
            pruned_snps.add(line.strip())

    # Filter the VCF file based on prune.in
    with gzip.open(input_vcf_gz, 'rt') as vcf_in, open(output_vcf, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
            else:
                fields = line.strip().split("\t")
                site = f"{fields[0]}:{fields[1]}"  # Chromosome:Position format
                if site in pruned_snps:
                    vcf_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <input.vcf.gz> <prune.in> <output.vcf>")
        sys.exit(1)

    input_vcf_gz = sys.argv[1]
    prune_in_file = sys.argv[2]
    output_vcf = sys.argv[3]

    LD_prunning_filtering(input_vcf_gz, prune_in_file, output_vcf)

