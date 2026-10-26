#!/usr/bin/env bash
# Non-model-organism track: do the polisher results (and Medaka's bacterial methylation model)
# generalize beyond model organisms? Two isolates from the Autocycler paper (Wick et al. 2025),
# chosen for phylogenetic distance from the usual training species:
#   Linnocua   Listeria innocua AUSMDU00097349  (Bacillota)      ONT SRR32486273, Illumina SRR32484996
#   Prettgeri  Providencia rettgeri AUSMDU00026122 (Morganellaceae) ONT SRR32486258, Illumina SRR32458788
# ONT reads are R10.4.1, Dorado sup v5.0.0. Truth = the paper's same-isolate Trycycler references
# (figshare article 28910099). Caveat: those references were polished with Medaka + Polypolish +
# pyPOLCA, so a shared polisher blind spot cannot be excluded.
# Runs the full 8-combination matrix per sample, plus two bacterial-model medaka arms per sample.
# Usage: autocycler_nonmodel_benchmark.sh <root_folder> [threads]
set -euo pipefail
command -v hashbrown > /dev/null || { echo "activate the HashBrown conda env first" >&2; exit 1; }
ROOT=${1:?usage: autocycler_nonmodel_benchmark.sh <root_folder> [threads]}
THREADS=${2:-16}
HARNESS=$(dirname "$0")/../benchmark_polishers.sh
RAW=$ROOT/data_autocycler_raw
DATA=$ROOT/data_autocycler
RUNS=$ROOT/runs_autocycler
mkdir -p "$RAW" "$DATA/Linnocua" "$DATA/Prettgeri"

# Reads from ENA
for path in \
    SRR324/073/SRR32486273/SRR32486273_1.fastq.gz \
    SRR324/096/SRR32484996/SRR32484996_1.fastq.gz \
    SRR324/096/SRR32484996/SRR32484996_2.fastq.gz \
    SRR324/058/SRR32486258/SRR32486258_1.fastq.gz \
    SRR324/088/SRR32458788/SRR32458788_1.fastq.gz \
    SRR324/088/SRR32458788/SRR32458788_2.fastq.gz ; do
    [ -s "$RAW/$(basename "$path")" ] || wget -q -c -P "$RAW" "https://ftp.sra.ebi.ac.uk/vol1/fastq/$path"
done
# References from the paper's figshare archive (contains all assemblies; we keep the references)
if [ ! -s "$RAW/assemblies.tar.xz" ]; then
    wget -q -c -O "$RAW/assemblies.tar.xz" https://ndownloader.figshare.com/files/56978396
fi
if [ ! -s "$DATA/Linnocua/reference.fasta" ]; then
    mkdir -p "$RAW/refs"
    tar -xJf "$RAW/assemblies.tar.xz" -C "$RAW/refs" --wildcards '*reference*.fasta'
    # one combined reference.fasta per species if present, otherwise the per-replicon parts
    pick_ref() { # species_glob -> stdout
        local combined
        combined=$(find "$RAW/refs" -path "*$1*" -name 'reference.fasta' | head -1)
        if [ -n "$combined" ]; then cat "$combined"
        else find "$RAW/refs" -path "*$1*" -name 'reference_*.fasta' | sort | xargs cat; fi
    }
    pick_ref Listeria > "$DATA/Linnocua/reference.fasta"
    pick_ref Providencia > "$DATA/Prettgeri/reference.fasta"
fi

# Stage: ONT filtered to >=1 kb best ~100x (genomes: 2.97 Mb and 4.47 Mb); Illumina subsampled
# to ~100x (raw depths ~270x and ~124x)
stage() { # sample ont_acc r_acc target_bases ill_frac
    local s=$1 ont=$2 ill=$3 tb=$4 frac=$5
    [ -s "$DATA/$s/$s.fastq.gz" ] || filtlong --min_length 1000 --keep_percent 95 \
        --target_bases "$tb" "$RAW/${ont}_1.fastq.gz" 2> /dev/null | gzip > "$DATA/$s/$s.fastq.gz"
    for mate in 1 2; do
        [ -s "$DATA/$s/${s}_R$mate.fastq.gz" ] || seqtk sample -s42 "$RAW/${ill}_$mate.fastq.gz" "$frac" \
            | gzip > "$DATA/$s/${s}_R$mate.fastq.gz"
    done
}
stage Linnocua  SRR32486273 SRR32484996 297000000 0.371
stage Prettgeri SRR32486258 SRR32458788 447000000 0.803

# 1. Full matrix (assembles each sample once, resumes 8 polisher combinations)
bash "$HARNESS" "$DATA" "$RUNS" "$THREADS"

# 2. Bacterial-model medaka arms per sample, resumed from each assembled base
field() { awk -F'\t' -v k="$2" '$1==k{print $2;f=1} END{if(!f)print "NA"}' "$1"; }
SUMMARY=$RUNS/bacterial_results.tsv
printf 'Sample\tLongPolisher\tShortPolisher\tMismatchesPer100kbp\tIndelsPer100kbp\tMisassemblies\tNGA50\tTotalLength\tContigs\n' > "$SUMMARY"
for s in Linnocua Prettgeri; do
    OUT=$RUNS/$s
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
        asm="$run/6_polished_short/$s.fasta"
        [ -s "$asm" ] || asm="$run/4_polished_long/$s.fasta"
        quast.py -r "$DATA/$s/reference.fasta" -o "$run/quast_vs_truth" -t "$THREADS" --no-icarus "$asm" \
            > /dev/null 2>&1
        r="$run/quast_vs_truth/report.tsv"
        printf '%s\tmedaka_bacterial\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$sp" \
            "$(field "$r" '# mismatches per 100 kbp')" "$(field "$r" '# indels per 100 kbp')" \
            "$(field "$r" '# misassemblies')" "$(field "$r" 'NGA50')" \
            "$(field "$r" 'Total length')" "$(field "$r" '# contigs')" >> "$SUMMARY"
    done
done
echo "Non-model benchmark done"
column -t "$RUNS/benchmark_results.tsv"
column -t "$SUMMARY"
