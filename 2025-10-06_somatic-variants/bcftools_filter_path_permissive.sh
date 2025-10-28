#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./bcftools_filter_path_permissive.sh -i input.vcf.gz -g genes.txt [-o {input}.filtered.vcf.gz]
#
# Filters:
#   - IMPACT == HIGH or MODERATE
#   - SYMBOL in genes.txt
#   - PolyPhen: possibly_damaging or probably_damaging (if present)
#   - SIFT: deleterious or deleterious_low_confidence (if present)
# Requires:
#   - bcftools
#   - split-vep plugin

# This is where the plugin will be in a conda environment.
export BCFTOOLS_PLUGINS=$CONDA_PREFIX/libexec/bcftools

VCF=""
GENES=""
OUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            VCF="$2"
            shift 2
            ;;
        -g|--genes)
            GENES="$2"
            shift 2
            ;;
        -o|--output)
            OUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 -i input.vcf.gz -g genes.txt [-o output.vcf.gz]"
            exit 1
            ;;
    esac
done

if [[ -z "$VCF" || -z "$GENES" ]]; then
    echo "Error: must specify both -i|--input and -g|--genes"
    exit 1
fi

# Default output if not specified
if [[ -z "$OUT" ]]; then
    base="$(basename "$VCF")"
    base="${base%.vcf.gz}"
    base="${base%.vcf.bgz}"
    base="${base%.vcf}"
    OUT="${base}.filtered.vcf.gz"
fi

echo "Filtering $VCF using gene list $GENES ..."
echo "Output: $OUT"

targets=$(mktemp)
bcftools +split-vep "$VCF" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%IMPACT\t%PolyPhen\t%SIFT\n' -d | \
awk -v genes="$GENES" '
BEGIN {
    IGNORECASE = 1
    while ((getline g < genes) > 0) gene[g] = 1
}
gene[$5] &&
($6 == "HIGH" || $6 == "MODERATE") &&
(($7 == "" || $7 ~ /damaging/)) &&
(($8 == "" || $8 ~ /deleterious/)) {
    print $1 "\t" $2
}' | sort -u > "$targets"

if [[ ! -s "$targets" ]]; then
    echo "No variants passed filters. Writing empty vcf."
    bcftools view -h -Oz -o "$OUT" "$VCF"
else
    echo "Writing likely pathogenic variants..."
    bcftools view -T "$targets" -Oz -o "$OUT" "$VCF"
fi

bcftools index -t "$OUT"
rm -f "$targets"
echo "Done. Filtered VCF written to: $OUT"
