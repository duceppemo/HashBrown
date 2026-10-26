#!/usr/bin/env bash
# Simulate the AFsim sample: Nanopore (badread) + Illumina (ART) reads from the finished
# M. bovis AF2122/97 genome (LT708304.1). The simulation source is the exact truth.
# badread and art are not part of the HashBrown environment (both on bioconda).
# Usage: afsim_make_data.sh <root_folder>   -> writes <root_folder>/data_afsim/AFsim/
set -euo pipefail
ROOT=${1:?usage: afsim_make_data.sh <root_folder>}
DIR=$ROOT/data_afsim/AFsim
mkdir -p "$DIR"
cd "$DIR"

[ -s reference.fasta ] || curl -s \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LT708304.1&rettype=fasta&retmode=text' \
    > reference.fasta

# ~100x Nanopore at Q20 (2026-era quality), nanopore2023 error/qscore models
badread simulate --reference reference.fasta --quantity 100x --length 15000,13000 \
    --identity 20,4 --error_model nanopore2023 --qscore_model nanopore2023 --seed 42 \
    2> badread.log | gzip > AFsim.fastq.gz

# ~60x Illumina HiSeq 2500 150 bp pairs
art_illumina -ss HS25 -i reference.fasta -p -l 150 -f 60 -m 400 -s 50 -rs 42 --noALN -o AFsim_R \
    > art.log 2>&1
gzip -f AFsim_R1.fq AFsim_R2.fq
mv AFsim_R1.fq.gz AFsim_R1.fastq.gz
mv AFsim_R2.fq.gz AFsim_R2.fastq.gz
echo "AFsim data ready: $DIR (run: bash validation/benchmark_polishers.sh $ROOT/data_afsim $ROOT/runs_afsim)"
