#!/usr/bin/env bash
# E. coli K-12 MG1655 (DSM 18039) against U00096.3: real reads on all three platforms.
# 1. ONT matrix: the full 8 polisher combinations (real R10.4.1 reads, default medaka model).
# 2. Two extra arms with ONT's bacterial methylation medaka model.
# 3. HiFi-base matrix: flye --pacbio-hifi base injected, same combinations — does polishing
#    damage an already excellent assembly?
# Truth caveat: the reads are from DSM 18039, the truth is U00096.3 — a handful of real lab-stock
# differences (~2 mismatches + 6 indels) form an error floor no polisher can cross.
# Needs sra-tools for the Illumina run (no ENA fastq mirror) and seqtk for subsampling.
# Usage: mg1655_benchmark.sh <root_folder> [threads]
set -euo pipefail
command -v hashbrown > /dev/null || { echo "activate the HashBrown conda env first" >&2; exit 1; }
ROOT=${1:?usage: mg1655_benchmark.sh <root_folder> [threads]}
THREADS=${2:-16}
HARNESS=$(dirname "$0")/../benchmark_polishers.sh
RAW=$ROOT/data_mg1655_raw
M=$ROOT/data_mg1655/MG1655
P=$ROOT/data_mg1655pb/MG1655pb
mkdir -p "$RAW" "$M" "$P"

# Truth: U00096.3
if [ ! -s "$M/reference.fasta" ]; then
    curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=fasta&retmode=text' \
        > "$M/reference.fasta"
fi
cp -f "$M/reference.fasta" "$P/reference.fasta"

# Reads: ONT R10.4.1 (ERR17025883, MicroBench) and HiFi (SRR10971019) from ENA,
# Illumina (SRR41214148) via sra-tools (ENA mirrors no fastq for it)
[ -s "$RAW/ERR17025883.fastq.gz" ] || wget -q -c -P "$RAW" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR170/083/ERR17025883/ERR17025883.fastq.gz"
[ -s "$RAW/SRR10971019_subreads.fastq.gz" ] || wget -q -c -P "$RAW" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/019/SRR10971019/SRR10971019_subreads.fastq.gz"
if [ ! -s "$RAW/SRR41214148_1.fastq" ]; then
    prefetch SRR41214148 --max-size 30G -O "$RAW/sra_tmp"
    fasterq-dump "$RAW/sra_tmp/SRR41214148" --threads "$THREADS" -O "$RAW"
fi

# Stage: long reads filtered to >=1 kb best ~100x; Illumina subsampled to ~60x (930k pairs)
[ -s "$M/MG1655.fastq.gz" ] || filtlong --min_length 1000 --keep_percent 95 \
    --target_bases 464200000 "$RAW/ERR17025883.fastq.gz" 2> /dev/null | gzip > "$M/MG1655.fastq.gz"
[ -s "$P/MG1655pb.fastq.gz" ] || filtlong --min_length 1000 --keep_percent 95 \
    --target_bases 464200000 "$RAW/SRR10971019_subreads.fastq.gz" 2> /dev/null | gzip > "$P/MG1655pb.fastq.gz"
for mate in 1 2; do
    [ -s "$M/MG1655_R$mate.fastq.gz" ] || seqtk sample -s42 "$RAW/SRR41214148_$mate.fastq" 930000 \
        | gzip > "$M/MG1655_R$mate.fastq.gz"
    cp -f "$M/MG1655_R$mate.fastq.gz" "$P/MG1655pb_R$mate.fastq.gz"
done

# 1. ONT matrix
bash "$HARNESS" "$ROOT/data_mg1655" "$ROOT/runs_mg1655" "$THREADS"

# 2. Bacterial-model arms, resumed from the assembled base
REF=$M/reference.fasta
OUT=$ROOT/runs_mg1655/MG1655
field() { awk -F'\t' -v k="$2" '$1==k{print $2;f=1} END{if(!f)print "NA"}' "$1"; }
SUMMARY=$ROOT/runs_mg1655/bacterial_results.tsv
printf 'Sample\tLongPolisher\tShortPolisher\tMismatchesPer100kbp\tIndelsPer100kbp\tMisassemblies\tNGA50\tTotalLength\tContigs\n' > "$SUMMARY"
for sp in none polypolish+pypolca; do
    run=$OUT/bacterial__${sp}
    if [ ! -d "$run" ]; then
        cp -r "$OUT/base" "$run"
        rm -rf "$run/done_polishing_long" "$run/done_trimming_short" "$run/done_polishing_short" \
               "$run/4_polished_long" "$run/5a_trimmed_short" "$run/5b_trimmed_short_report" \
               "$run/6_polished_short" "$run/7_assembly_qc"
    fi
    hashbrown -l "$OUT/reads_long" -s "$OUT/reads_short" -o "$run" \
        --long-polisher medaka --model r1041_e82_400bps_bacterial_methylation \
        --short-polisher "$sp" -t "$THREADS"
    asm="$run/6_polished_short/MG1655.fasta"
    [ -s "$asm" ] || asm="$run/4_polished_long/MG1655.fasta"
    quast.py -r "$REF" -o "$run/quast_vs_truth" -t "$THREADS" --no-icarus "$asm" > /dev/null 2>&1
    r="$run/quast_vs_truth/report.tsv"
    printf 'MG1655\tmedaka_bacterial\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sp" \
        "$(field "$r" '# mismatches per 100 kbp')" "$(field "$r" '# indels per 100 kbp')" \
        "$(field "$r" '# misassemblies')" "$(field "$r" 'NGA50')" \
        "$(field "$r" 'Total length')" "$(field "$r" '# contigs')" >> "$SUMMARY"
done
column -t "$SUMMARY"

# 3. HiFi-base matrix
base=$ROOT/runs_mg1655pb/MG1655pb/base
if [ ! -s "$base/3_assembled_long/MG1655pb.fasta" ]; then
    flye --pacbio-hifi "$P/MG1655pb.fastq.gz" --out-dir "$ROOT/flye_hifi" --threads "$THREADS"
    mkdir -p "$base/3_assembled_long"
    cp "$ROOT/flye_hifi/assembly.fasta" "$base/3_assembled_long/MG1655pb.fasta"
    touch "$base/done_assembling_long"
fi
bash "$HARNESS" "$ROOT/data_mg1655pb" "$ROOT/runs_mg1655pb" "$THREADS"
echo "MG1655 benchmarks done"
