#!/usr/bin/env bash
# Medaka-vs-none on real R10.4.1 sup v4.3.0 Nanopore reads (M. bovis field strain MBWGS426),
# scored against the AF2122/97 genome. Field strain vs reference: strain divergence dominates the
# absolute numbers, so compare the arms with each other, not with zero.
# The reads are in-house (not public): point READS at your own R10.4.1 bacterial fastq to rerun.
# Usage: mbwgs426_medaka_models.sh <root_folder> <reads.fastq.gz> [threads]
set -euo pipefail
command -v hashbrown > /dev/null || { echo "activate the HashBrown conda env first" >&2; exit 1; }
ROOT=${1:?usage: mbwgs426_medaka_models.sh <root_folder> <reads.fastq.gz> [threads]}
READS=${2:?path to the Nanopore reads of the sample}
THREADS=${3:-16}
OUT=$ROOT/mbwgs426
SAMPLE=$(basename "$READS" | sed 's/\.fastq\.gz$//; s/\.fq\.gz$//')

REF=$OUT/reference.fasta
mkdir -p "$OUT/reads_long"
[ -s "$REF" ] || curl -s \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LT708304.1&rettype=fasta&retmode=text' \
    > "$REF"
ln -sf "$(realpath "$READS")" "$OUT/reads_long/$(basename "$READS")"

# Assemble once, no polishing
hashbrown -l "$OUT/reads_long" -o "$OUT/none" --long-polisher none -t "$THREADS" -p 1

declare -A MODELS=(
  [medaka_sup_v4.3.0]=r1041_e82_400bps_sup_v4.3.0
  [medaka_bacterial]=r1041_e82_400bps_bacterial_methylation
)
for arm in "${!MODELS[@]}"; do
    if [ ! -d "$OUT/$arm" ]; then
        cp -r "$OUT/none" "$OUT/$arm"
        rm -rf "$OUT/$arm/4_polished_long" "$OUT/$arm/done_polishing_long" "$OUT/$arm/7_assembly_qc"
    fi
    hashbrown -l "$OUT/reads_long" -o "$OUT/$arm" --model "${MODELS[$arm]}" -t "$THREADS" -p 1
done

SUMMARY=$OUT/mbwgs426_results.tsv
printf 'Arm\tMismatchesPer100kbp\tIndelsPer100kbp\tMisassemblies\tNGA50\tTotalLength\tContigs\n' > "$SUMMARY"
field() { awk -F'\t' -v k="$2" '$1==k{print $2;f=1} END{if(!f)print "NA"}' "$1"; }
for arm in none medaka_sup_v4.3.0 medaka_bacterial; do
    asm="$OUT/$arm/4_polished_long/$SAMPLE.fasta"
    [ -s "$asm" ] || asm="$OUT/$arm/3_assembled_long/$SAMPLE.fasta"
    quast.py -r "$REF" -o "$OUT/$arm/quast_vs_truth" -t "$THREADS" --no-icarus "$asm" > /dev/null 2>&1
    r="$OUT/$arm/quast_vs_truth/report.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$arm" \
        "$(field "$r" '# mismatches per 100 kbp')" "$(field "$r" '# indels per 100 kbp')" \
        "$(field "$r" '# misassemblies')" "$(field "$r" 'NGA50')" \
        "$(field "$r" 'Total length')" "$(field "$r" '# contigs')" >> "$SUMMARY"
done
echo "MBWGS426 benchmark done"
column -t "$SUMMARY"
