"""Stubs for the external programs, so the pipeline tests run in seconds without conda."""

import gzip
import os
from pathlib import Path

import pytest

FASTA = '>contig_1\nACGTACGTACGT\n'
SAM = 'read1\t0\tcontig_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'

# Each stub mimics just enough of the real program's command line for hashbrown's calls.
STUBS = {
    'porechop': '''
    if [ "$1" = "--version" ]; then echo 0.2.4; exit 0; fi
    in=""; out=""; prev=""
    for a in "$@"; do [ "$prev" = "-i" ] && in=$a; [ "$prev" = "-o" ] && out=$a; prev=$a; done
    cp "$in" "$out"
    ''',
    'filtlong': '''
    if [ "$1" = "--version" ]; then echo "Filtlong v0.2.1"; exit 0; fi
    for a in "$@"; do in=$a; done
    case "$in" in *.gz) gzip -dc "$in";; *) cat "$in";; esac
    ''',
    'flye': '''
    if [ "$1" = "--version" ]; then echo 2.9.2; exit 0; fi
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--out-dir" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    printf '%s' ">contig_1\\nACGTACGTACGT\\n" > "$out/assembly.fasta"
    printf 'S\\tcontig_1\\tACGT\\n' > "$out/assembly_graph.gfa"
    {
      echo "Total read length: 120000"
      echo "Reads N50/N90: 8000/3000"
      printf 'Total length:\\t12000\\n'
      printf 'Fragments:\\t1\\n'
      printf 'Mean coverage:\\t50\\n'
    } > "$out/flye.log"
    ''',
    'shasta': '''
    if [ "$1" = "--version" ]; then echo "Shasta Release 0.11.1"; exit 0; fi
    out=""; prev=""; cleanup=0
    for a in "$@"; do
      [ "$prev" = "--assemblyDirectory" ] && out=$a
      [ "$prev" = "--command" ] && cleanup=1
      prev=$a
    done
    [ "$cleanup" = 1 ] && exit 0
    mkdir -p "$out"
    printf '%s' ">contig_1\\nACGTACGTACGT\\n" > "$out/Assembly.fasta"
    printf 'S\\tcontig_1\\tACGT\\n' > "$out/Assembly.gfa"
    echo "Total number of reads is 1000."
    echo "Total number of raw bases is 9000000."
    echo "N50 for read length is 8000 bases."
    echo "Total length of assembled sequence is 12000"
    echo "The assembly graph has 1 assembled segments"
    ''',
    'medaka_consensus': '''
    if [ "$1" = "-h" ]; then echo "    Choices: r1041_e82_400bps_sup_g615 model_a model_b"; exit 0; fi
    asm=""; out=""; prev=""
    for a in "$@"; do [ "$prev" = "-d" ] && asm=$a; [ "$prev" = "-o" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    cp "$asm" "$out/consensus.fasta"
    ''',
    'medaka': 'echo "medaka 2.0.1"',
    'fastp': '''
    in1=""; in2=""; out1=""; out2=""; html=""; json=""
    while [ $# -gt 0 ]; do
      case "$1" in
        --version) echo "fastp 0.23.4"; exit 0;;
        --in1) in1=$2; shift;; --in2) in2=$2; shift;;
        --out1) out1=$2; shift;; --out2) out2=$2; shift;;
        --html) html=$2; shift;; --json) json=$2; shift;;
      esac; shift
    done
    [ -n "$out1" ] && cp "$in1" "$out1"
    [ -n "$out2" ] && [ -n "$in2" ] && cp "$in2" "$out2"
    [ -n "$html" ] && : > "$html"
    [ -n "$json" ] && : > "$json"
    exit 0
    ''',
    'bwa': '''
    if [ "$1" = "index" ]; then for e in amb ann bwt pac sa; do : > "$2.$e"; done
    elif [ "$1" = "mem" ]; then printf '%s' "''' + SAM + '''"
    fi
    ''',
    'polypolish': '''
    if [ "$1" = "--version" ]; then echo "Polypolish v0.6.0"; exit 0; fi
    if [ "$1" = "polish" ]; then cat "$2"; echo "polypolish log" >&2; fi
    ''',
    'pypolca': '''
    asm=""; out=""; prefix=""
    while [ $# -gt 0 ]; do
      case "$1" in
        --version) echo "pypolca 0.3.1"; exit 0;;
        --assembly) asm=$2; shift;; --output) out=$2; shift;; --prefix) prefix=$2; shift;;
      esac; shift
    done
    mkdir -p "$out"
    cp "$asm" "$out/${prefix}_corrected.fasta"
    echo report > "$out/$prefix.report"
    echo vcf > "$out/$prefix.vcf"
    echo "pypolca done"
    ''',
    'skesa': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--contigs_out" ] && out=$a; prev=$a; done
    [ -n "$out" ] && printf '%s' ">contig_1\\nACGTACGTACGT\\n" > "$out"
    exit 0
    ''',
    'spades.py': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    printf '%s' ">NODE_1\\nACGTACGTACGT\\n" > "$out/scaffolds.fasta"
    printf 'S\\t1\\tACGT\\n' > "$out/assembly_graph_with_scaffolds.gfa"
    ''',
    'gfa_connector': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--gfa" ] && out=$a; prev=$a; done
    [ -n "$out" ] && printf 'S\\tcontig_1\\tACGT\\n' > "$out"
    exit 0
    ''',
    'seqtk': '''
    if [ "$1" = "seq" ]; then
      shift
      if [ "$1" = "-a" ]; then printf '%s' ">read1\\nACGTACGTACGT\\n"; else cat "$1"; fi
    fi
    exit 0
    ''',
    'minimap2': '''
    if [ "$1" = "--version" ]; then echo 2.26-r1175; exit 0; fi
    printf '%s' "''' + SAM + '''"
    ''',
    'samtools': '''
    case "$1" in
      view|sort|fixmate) cat;;
      markdup) for a in "$@"; do out=$a; done; cat > "$out";;
      index) : > "$2.bai";;
      depth) printf 'contig_1\\t1\\t10\\ncontig_1\\t2\\t10\\ncontig_1\\t3\\t10\\n';;
      version) echo "samtools 1.17";;
    esac
    exit 0
    ''',
    'qualimap': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-outdir" ] && out=$a; prev=$a; done
    if [ -z "$out" ]; then echo "QualiMap v.2.2.2"; exit 0; fi
    mkdir -p "$out"
    echo report > "$out/qualimapReport.html"
    ''',
    'quast.py': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--output-dir" ] && out=$a; prev=$a; done
    if [ -z "$out" ]; then echo "QUAST v5.2.0"; exit 0; fi
    mkdir -p "$out"
    printf 'Assembly\\tsample\\n' > "$out/report.tsv"
    ''',
    'nucmer': '''
    if [ "$1" = "--version" ]; then echo 4.0.0; exit 0; fi
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-p" ] && out=$a; prev=$a; done
    [ -n "$out" ] && : > "$out.delta"
    exit 0
    ''',
    'mummerplot': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-p" ] && out=$a; prev=$a; done
    [ -n "$out" ] && : > "$out.ps"
    exit 0
    ''',
    'gnuplot': 'exit 0',
    'Bandage': '''
    if [ "$1" = "--version" ]; then echo "Version: 0.8.1"; exit 0; fi
    [ "$1" = "image" ] && : > "$3"
    exit 0
    ''',
}


@pytest.fixture
def stub_bin(tmp_path_factory, monkeypatch):
    """A folder of fake external programs, put first on the PATH."""
    bin_dir = tmp_path_factory.mktemp('stub_bin')
    for name, body in STUBS.items():
        stub = bin_dir / name
        stub.write_text('#!/usr/bin/env bash\n' + '\n'.join(line.strip() for line in body.strip().splitlines())
                        + '\n')
        stub.chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')
    return bin_dir


def write_fastq_gz(path: Path, n_reads: int = 3) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt') as f:
        for i in range(n_reads):
            f.write(f'@read{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n')
    return path


@pytest.fixture
def fastq_gz():
    return write_fastq_gz
