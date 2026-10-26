#!/usr/bin/env bash
# Polisher matrix on real M. bovis AF2122/97 reads (PacBio RS CLR + NextSeq pair) against the
# strain's own finished genome (LT708304.1): real reads, exact truth. The PacBio reads are CLR, so
# the base assembly is made with flye --pacbio-raw outside hashbrown and injected into the harness
# state (hashbrown's resume skips the assembly step when its flag and output are present).
# Caveat for the report: Medaka is Nanopore-trained; its arms on PacBio input are indicative only.
# Needs sra-tools (prefetch/fasterq-dump) for the PacBio run; Illumina comes from ENA.
# Usage: af2122_benchmark.sh <root_folder> [threads]
set -euo pipefail
command -v hashbrown > /dev/null || { echo "activate the HashBrown conda env first" >&2; exit 1; }
ROOT=${1:?usage: af2122_benchmark.sh <root_folder> [threads]}
THREADS=${2:-16}
HARNESS=$(dirname "$0")/../benchmark_polishers.sh
DATA=$ROOT/data_af2122
RAW=$ROOT/data_af2122_raw
RUNS=$ROOT/runs_af2122
mkdir -p "$DATA/AF2122" "$RAW"

# Truth: the strain's own finished genome
[ -s "$DATA/AF2122/reference.fasta" ] || curl -s \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LT708304.1&rettype=fasta&retmode=text' \
    > "$DATA/AF2122/reference.fasta"

# PacBio RS CLR subreads (ERR1748639; original submission is PacBio-native, hence sra-tools)
if [ ! -s "$RAW/ERR1748639.fastq.gz" ]; then
    prefetch ERR1748639 --max-size 50G -O "$RAW/sra_tmp"
    fasterq-dump "$RAW/sra_tmp/ERR1748639" --threads "$THREADS" -O "$RAW"
    gzip "$RAW/ERR1748639.fastq"
fi
# NextSeq pair (ERR4586795) from ENA
for mate in 1 2; do
    [ -s "$RAW/ERR4586795_$mate.fastq.gz" ] || wget -q -c -P "$RAW" \
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/005/ERR4586795/ERR4586795_$mate.fastq.gz"
done

# Stage: CLR filtered to >=1 kb, best ~100x; Illumina linked as the sample's pair
if [ ! -s "$DATA/AF2122/AF2122.fastq.gz" ]; then
    filtlong --min_length 1000 --keep_percent 95 --target_bases 435000000 \
        "$RAW/ERR1748639.fastq.gz" 2> /dev/null | gzip > "$DATA/AF2122/AF2122.fastq.gz"
fi
ln -sf "$RAW/ERR4586795_1.fastq.gz" "$DATA/AF2122/AF2122_R1.fastq.gz"
ln -sf "$RAW/ERR4586795_2.fastq.gz" "$DATA/AF2122/AF2122_R2.fastq.gz"

# Assemble the CLR reads with the right Flye mode and inject the result as the base run
base=$RUNS/AF2122/base
if [ ! -s "$base/3_assembled_long/AF2122.fasta" ]; then
    flye --pacbio-raw "$DATA/AF2122/AF2122.fastq.gz" --out-dir "$ROOT/flye_pacbio" \
         --threads "$THREADS" --iterations 3
    mkdir -p "$base/3_assembled_long"
    cp "$ROOT/flye_pacbio/assembly.fasta" "$base/3_assembled_long/AF2122.fasta"
    touch "$base/done_assembling_long"
fi

bash "$HARNESS" "$DATA" "$RUNS" "$THREADS"
