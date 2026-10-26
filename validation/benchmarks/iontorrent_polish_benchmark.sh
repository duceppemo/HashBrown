#!/usr/bin/env bash
# Does polishing a Nanopore assembly with Ion Torrent reads work, or does it fail on exactly the
# homopolymer indels Illumina fixes (both platforms share that error mode)?
#
# Two isolates, each sequenced on THREE platforms from one DNA submission:
#   Sgoete   Scandinavium goeteborgense CCUG 66741 (type strain)
#            ONT SRR11858646 (R9.4.1 MinION), Ion PGM SRR8561413 (SE), NovaSeq SRR11858647 (PE)
#   VibED400 Vibrio sp. ED-400
#            ONT SRR16167661 (MinION), Ion PGM SRR16167662 (SE), MiSeq SRR16167664 (PE)
#
# Design:
#   truth  = Flye on the BEST ~100x ONT + Polypolish + pyPOLCA with the FULL Illumina run
#            (no Medaka: the R9 basecaller versions are unknown, and the question here is only
#            about the short-read platform). Caveat: the truth is Illumina-polished, which favors
#            the Illumina arm; the Ion-vs-none comparison and the homopolymer classification of
#            the remaining errors are the meaningful readouts.
#   base   = Flye on a RANDOM ~50x ONT subsample (worse than the truth's base on purpose:
#            polishing needs something to fix).
#   arms   = base + nothing | + Ion (SE) | + Illumina subsampled to the Ion arm's base count.
#            Both polish arms run Polypolish then pyPOLCA --careful.
# Scoring: QUAST vs truth, plus every remaining indel classified by homopolymer context
#          (a run of >=3 identical bases around the indel in the truth).
# Usage: iontorrent_polish_benchmark.sh <root_folder> [threads]
set -euo pipefail
command -v hashbrown > /dev/null || { echo "activate the HashBrown conda env first" >&2; exit 1; }
ROOT=${1:?usage: iontorrent_polish_benchmark.sh <root_folder> [threads]}
THREADS=${2:-16}
RAW=$ROOT/data_iontorrent_raw
RUNS=$ROOT/runs_iontorrent
mkdir -p "$RAW" "$RUNS"

ena() { [ -s "$RAW/$(basename "$1")" ] || wget -q -c -P "$RAW" "https://ftp.sra.ebi.ac.uk/vol1/fastq/$1"; }
ena SRR118/046/SRR11858646/SRR11858646_1.fastq.gz   # Sgoete ONT
ena SRR856/003/SRR8561413/SRR8561413.fastq.gz       # Sgoete Ion (SE)
ena SRR118/047/SRR11858647/SRR11858647_1.fastq.gz   # Sgoete NovaSeq R1
ena SRR118/047/SRR11858647/SRR11858647_2.fastq.gz   # Sgoete NovaSeq R2
ena SRR161/061/SRR16167661/SRR16167661_1.fastq.gz   # VibED400 ONT
ena SRR161/062/SRR16167662/SRR16167662.fastq.gz     # VibED400 Ion (SE)
ena SRR161/064/SRR16167664/SRR16167664_1.fastq.gz   # VibED400 MiSeq R1
ena SRR161/064/SRR16167664/SRR16167664_2.fastq.gz   # VibED400 MiSeq R2

bases_of() { seqkit stats -T "$1" | awk 'NR==2{print $5}'; }

field() { awk -F'\t' -v k="$2" '$1==k{print $2;f=1} END{if(!f)print "NA"}' "$1"; }

SUMMARY=$RUNS/iontorrent_results.tsv
printf 'Sample\tArm\tMismatchesPer100kbp\tIndelsPer100kbp\tIndels\tHomopolymerIndels\tMisassemblies\n' \
    > "$SUMMARY"

run_sample() { # name ont_fastq ion_fastq ill_r1 ill_r2
    local name=$1 ont=$RAW/$2 ion=$RAW/$3 r1=$RAW/$4 r2=$RAW/$5
    local work=$RUNS/$name
    mkdir -p "$work"

    # Truth: best ~100x ONT + full Illumina polish
    if [ ! -s "$work/truth.fasta" ]; then
        mkdir -p "$work/truth_long" "$work/truth_short"
        [ -s "$work/truth_long/$name.fastq.gz" ] || filtlong --min_length 1000 --keep_percent 95 \
            --target_bases 500000000 "$ont" 2> /dev/null | gzip > "$work/truth_long/$name.fastq.gz"
        ln -sf "$r1" "$work/truth_short/${name}_R1.fastq.gz"
        ln -sf "$r2" "$work/truth_short/${name}_R2.fastq.gz"
        hashbrown -l "$work/truth_long" -s "$work/truth_short" -o "$work/truth_run" \
            --long-polisher none -t "$THREADS" -p 1 --skip-report
        cp "$work/truth_run/6_polished_short/$name.fasta" "$work/truth.fasta"
    fi

    # Arms base: random ~50x ONT subsample, assembled once without polishing
    if [ ! -s "$work/arms/none/3_assembled_long/$name.fasta" ]; then
        mkdir -p "$work/arm_long"
        if [ ! -s "$work/arm_long/$name.fastq.gz" ]; then
            local ont_bases frac
            ont_bases=$(bases_of "$ont")
            frac=$(python3 -c "print(min(1.0, 250000000 / $ont_bases))")
            seqtk sample -s7 "$ont" "$frac" | gzip > "$work/arm_long/$name.fastq.gz"
        fi
        mkdir -p "$work/arms"
        hashbrown -l "$work/arm_long" -o "$work/arms/none" --long-polisher none \
            --short-polisher none -t "$THREADS" -p 1 --skip-report
    fi

    # Illumina subsample matched to the Ion base count (fair fight)
    if [ ! -s "$work/ill_matched/${name}_R1.fastq.gz" ]; then
        mkdir -p "$work/ill_matched"
        local ion_bases ill_bases frac
        ion_bases=$(bases_of "$ion")
        ill_bases=$(($(bases_of "$r1") + $(bases_of "$r2")))
        frac=$(python3 -c "print(min(1.0, $ion_bases / $ill_bases))")
        seqtk sample -s7 "$r1" "$frac" | gzip > "$work/ill_matched/${name}_R1.fastq.gz"
        seqtk sample -s7 "$r2" "$frac" | gzip > "$work/ill_matched/${name}_R2.fastq.gz"
    fi
    mkdir -p "$work/ion_reads"
    ln -sf "$ion" "$work/ion_reads/$name.fastq.gz"

    # Polish arms, resumed from the shared base
    local arm
    for arm in ion illumina; do
        if [ ! -d "$work/arms/$arm" ]; then
            cp -r "$work/arms/none" "$work/arms/$arm"
            rm -rf "$work/arms/$arm/done_polishing_short" "$work/arms/$arm/6_polished_short" \
                   "$work/arms/$arm/7_assembly_qc" "$work/arms/$arm/8_report" \
                   "$work/arms/$arm/done_reorienting"
        fi
        local short_dir=$work/ion_reads extra=(--ion-torrent)
        if [ "$arm" = illumina ]; then short_dir=$work/ill_matched; extra=(); fi
        hashbrown -l "$work/arm_long" -s "$short_dir" -o "$work/arms/$arm" \
            --long-polisher none --short-polisher polypolish+pypolca "${extra[@]}" \
            -t "$THREADS" -p 1 --skip-report
    done

    # Score every arm against the truth
    for arm in none ion illumina; do
        local asm=$work/arms/$arm/6_polished_short/$name.fasta
        [ -s "$asm" ] || asm=$work/arms/$arm/3_assembled_long/$name.fasta
        quast.py -r "$work/truth.fasta" -o "$work/arms/$arm/quast_vs_truth" -t "$THREADS" \
            --no-icarus "$asm" > /dev/null 2>&1
        nucmer -t "$THREADS" -p "$work/arms/$arm/vs_truth" "$work/truth.fasta" "$asm" 2> /dev/null
        show-snps -ClrT "$work/arms/$arm/vs_truth.delta" > "$work/arms/$arm/vs_truth.snps" || true
        local counts
        counts=$(python3 - "$work/truth.fasta" "$work/arms/$arm/vs_truth.snps" <<'PYEOF'
import sys
seqs, name = {}, None
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
seqs = {n: ''.join(p) for n, p in seqs.items()}
indels = homopolymer = 0
with open(sys.argv[2]) as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11 or not fields[0].isdigit():
            continue
        pos, ref_base, alt_base, contig = int(fields[0]), fields[1], fields[2], fields[10]
        if ref_base != '.' and alt_base != '.':
            continue  # mismatch, counted by QUAST
        indels += 1
        seq = seqs.get(contig, '')
        base = ref_base if ref_base != '.' else alt_base
        start = max(0, pos - 4)
        context = seq[start:pos + 4]
        if base * 3 in context:
            homopolymer += 1
print(f'{indels}\t{homopolymer}')
PYEOF
)
        local r=$work/arms/$arm/quast_vs_truth/report.tsv
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$arm" \
            "$(field "$r" '# mismatches per 100 kbp')" "$(field "$r" '# indels per 100 kbp')" \
            "$counts" "$(field "$r" '# misassemblies')" >> "$SUMMARY"
    done
}

run_sample Sgoete SRR11858646_1.fastq.gz SRR8561413.fastq.gz \
    SRR11858647_1.fastq.gz SRR11858647_2.fastq.gz
run_sample VibED400 SRR16167661_1.fastq.gz SRR16167662.fastq.gz \
    SRR16167664_1.fastq.gz SRR16167664_2.fastq.gz

echo "Ion Torrent polishing benchmark done"
column -t "$SUMMARY"
