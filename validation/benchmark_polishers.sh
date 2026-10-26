#!/usr/bin/env bash
# Benchmark every polisher combination against gold-standard references.
#
# For each sample, the long reads are assembled once, then each combination of
#   --long-polisher  {none, medaka}
#   --short-polisher {none, polypolish, pypolca, polypolish+pypolca}
# is applied by copying the assembled state and resuming, so the assembly is not recomputed.
# Each final assembly is compared with the gold standard by QUAST (-r), and the error rates
# (mismatches and indels per 100 kbp, misassemblies, NGA50) land in one table.
#
# Usage: bash validation/benchmark_polishers.sh <data_folder> <output_folder> [threads]
#
# <data_folder> holds one subfolder per sample:
#   <sample>/reference.fasta        the gold standard (a finished, trusted genome)
#   <sample>/<sample>.fastq.gz      the Nanopore reads (name must start with "<sample>_" or "<sample>.")
#   <sample>/<sample>_R1.fastq.gz   the Illumina reads
#   <sample>/<sample>_R2.fastq.gz
#
# Reads can be real (sequenced from the same culture as the reference) or simulated from the
# reference (e.g. badread for Nanopore, wgsim or art_illumina for Illumina), in which case the
# reference is the exact truth.

set -euo pipefail

DATA=${1:?data folder}
OUT=${2:?output folder}
THREADS=${3:-$(nproc)}

LONG_POLISHERS=(none medaka)
SHORT_POLISHERS=(none polypolish pypolca polypolish+pypolca)
# Extra hashbrown options for every run, e.g. HASHBROWN_EXTRA_ARGS="--model r1041_e82_400bps_sup_v4.3.0"
EXTRA_ARGS=${HASHBROWN_EXTRA_ARGS:-}

command -v hashbrown >/dev/null || { echo "hashbrown is not on the PATH" >&2; exit 1; }
command -v quast.py >/dev/null || { echo "quast.py is not on the PATH" >&2; exit 1; }

mkdir -p "$OUT"
SUMMARY="$OUT/benchmark_results.tsv"
printf 'Sample\tLongPolisher\tShortPolisher\tMismatchesPer100kbp\tIndelsPer100kbp\tMisassemblies\tNGA50\tTotalLength\tContigs\n' > "$SUMMARY"

final_assembly() {  # the most polished assembly of a run
    local run=$1 sample=$2
    for f in "$run/6_polished_short/$sample.fasta" \
             "$run/4_polished_long/$sample.fasta" \
             "$run/3_assembled_long/$sample.fasta"; do
        [ -s "$f" ] && { echo "$f"; return; }
    done
    return 1
}

quast_field() {  # value of one row of QUAST's report.tsv (second column)
    awk -F'\t' -v key="$2" '$1 == key { print $2; found = 1 } END { if (!found) print "NA" }' "$1"
}

for sample_dir in "$DATA"/*/; do
    sample=$(basename "$sample_dir")
    reference="$sample_dir/reference.fasta"
    [ -s "$reference" ] || { echo "No reference.fasta for $sample, skipping" >&2; continue; }

    # Stage the reads the way hashbrown expects them (one folder per technology). The link keeps
    # the data folder's file name even when it is itself a symlink to a file named differently.
    long_dir="$OUT/$sample/reads_long"; short_dir="$OUT/$sample/reads_short"
    mkdir -p "$long_dir" "$short_dir"
    for f in "$sample_dir/$sample"*.f*q.gz; do
        name=$(basename "$f")
        case "$name" in
            *_R1*|*_R2*) ln -sf "$(realpath "$f")" "$short_dir/$name" ;;
            *) ln -sf "$(realpath "$f")" "$long_dir/$name" ;;
        esac
    done

    # Assemble once (no polishing), then copy that state for each combination and resume
    base="$OUT/$sample/base"
    echo "=== $sample: assembling (base run) ==="
    hashbrown -l "$long_dir" -s "$short_dir" -o "$base" \
        --long-polisher none --short-polisher none --skip-report --no-trim-long --no-filter-long --no-trim-short -t "$THREADS" $EXTRA_ARGS

    for lp in "${LONG_POLISHERS[@]}"; do
        for sp in "${SHORT_POLISHERS[@]}"; do
            combo="${lp}__${sp}"
            run="$OUT/$sample/$combo"
            echo "=== $sample: $lp + $sp ==="
            if [ ! -d "$run" ]; then
                cp -r "$base" "$run"
                rm -rf "$run/8_report" "$run/done_polishing_long" "$run/done_trimming_short" "$run/done_polishing_short" \
                       "$run/4_polished_long" "$run/5a_trimmed_short" "$run/5b_trimmed_short_report" \
                       "$run/6_polished_short" "$run/7_assembly_qc"
            fi
            hashbrown -l "$long_dir" -s "$short_dir" -o "$run" \
                --long-polisher "$lp" --short-polisher "$sp" --skip-report --no-trim-long --no-filter-long --no-trim-short -t "$THREADS" $EXTRA_ARGS

            assembly=$(final_assembly "$run" "$sample") || { echo "No assembly for $sample/$combo" >&2; continue; }
            quast.py -r "$reference" -o "$run/quast_vs_truth" -t "$THREADS" --no-icarus "$assembly" \
                > /dev/null 2>&1 || true
            report="$run/quast_vs_truth/report.tsv"
            [ -s "$report" ] || { echo "QUAST failed for $sample/$combo" >&2; continue; }

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$sample" "$lp" "$sp" \
                "$(quast_field "$report" '# mismatches per 100 kbp')" \
                "$(quast_field "$report" '# indels per 100 kbp')" \
                "$(quast_field "$report" '# misassemblies')" \
                "$(quast_field "$report" 'NGA50')" \
                "$(quast_field "$report" 'Total length')" \
                "$(quast_field "$report" '# contigs')" >> "$SUMMARY"
        done
    done
done

echo
echo "Results: $SUMMARY (lower mismatches/indels per 100 kbp is better; the best combination wins)"
column -t "$SUMMARY"
