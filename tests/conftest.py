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
    'mosdepth': '''
    if [ "$1" = "--version" ]; then echo "mosdepth 0.3.11"; exit 0; fi
    prefix=""
    for a in "$@"; do prefix=$prev; prev=$a; done
    printf 'chrom\tlength\tbases\tmean\tmin\tmax\n' > "$prefix.mosdepth.summary.txt"
    printf 'contig_1\t12\t120\t10.00\t10\t10\n' >> "$prefix.mosdepth.summary.txt"
    printf 'total\t12\t120\t10.00\t10\t10\n' >> "$prefix.mosdepth.summary.txt"
    : > "$prefix.mosdepth.global.dist.txt"
    ''',
    'quast.py': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--output-dir" ] && out=$a; prev=$a; done
    if [ -z "$out" ]; then echo "QUAST v5.2.0"; exit 0; fi
    mkdir -p "$out"
    printf 'Assembly\\tsample\\n' > "$out/report.tsv"
    ''',
    'dnaapler': '''
    out=""; prefix=""; prev=""
    in=""
    for a in "$@"; do
      [ "$prev" = "-i" ] && in=$a; [ "$prev" = "-o" ] && out=$a; [ "$prev" = "-p" ] && prefix=$a
      prev=$a
    done
    mkdir -p "$out"
    cp "$in" "$out/${prefix}_reoriented.fasta"
    ''',
    'checkm2': '''
    [ "$1" = "database" ] && exit 0
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    printf 'Name\tCompleteness\tContamination\n' > "$out/quality_report.tsv"
    printf 'asm\t99.98\t0.12\n' >> "$out/quality_report.tsv"
    ''',
    'busco': '''
    out=""; path=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; [ "$prev" = "--out_path" ] && path=$a; prev=$a; done
    mkdir -p "$path/$out"
    {
      echo "# The lineage dataset is: bacteria_odb10 (Creation date: 2024-01-08)"
      echo "C:99.2%[S:99.0%,D:0.2%],F:0.3%,M:0.5%,n:124"
    } > "$path/$out/short_summary.specific.bacteria_odb10.txt"
    ''',
    'barrnap': '''
    if [ "$1" = "--version" ]; then echo "barrnap 0.9"; exit 0; fi
    printf 'contig_1\tbarrnap:0.9\trRNA\t1\t120\t0\t+\t.\tName=16S_rRNA;product=16S\n'
    printf 'contig_1\tbarrnap:0.9\trRNA\t200\t320\t0\t+\t.\tName=23S_rRNA;product=23S\n'
    printf 'contig_1\tbarrnap:0.9\trRNA\t400\t450\t0\t+\t.\tName=5S_rRNA;product=5S\n'
    ''',
    'tRNAscan-SE': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    {
      printf 'Sequence\ttRNA\tBounds\tBounds\ttRNA\tCodon\tIntron\tIntron\tInf\n'
      printf 'Name\t#\tBegin\tEnd\tType\t\t\t\tScore\n'
      printf -- '--------\t---\t----\t----\t----\t-----\t-----\t----\t-----\n'
      printf 'contig_1\t1\t100\t176\tAla\tAGC\t0\t0\t80.0\n'
      printf 'contig_1\t2\t300\t376\tGly\tGCC\t0\t0\t80.0\n'
    } > "$out"
    ''',
    'meryl': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "output" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    ''',
    'merqury.sh': '''
    printf 'asm\t100\t1000\t54.50\t0.00001\n' > "$3.qv"
    ''',
    'mash': '''
    if [ "$1" = "--version" ]; then echo 2.3; exit 0; fi
    printf 'refseq-NZ-GCF_000005845.2-Escherichia_coli.fna\tquery\t0.0010\t0\t995/1000\n'
    ''',
    'datasets': '''
    if [ "$1" = "summary" ]; then
      printf '{"reports":[{"organism":{"organism_name":"Escherichia coli"}}]}\n'
      exit 0
    fi
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--filename" ] && out=$a; prev=$a; done
    python3 -c "import sys,zipfile; z=zipfile.ZipFile(sys.argv[1],'w'); z.writestr('ncbi_dataset/data/GCF_000005845.2/ref.fna','>r\\nACGTACGT\\n'); z.close()" "$out"
    ''',
    'fastANI': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    printf 'q\tr\t99.99\t50\t50\n' > "$out"
    ''',
    'mlst': '''
    for a in "$@"; do in=$a; done
    printf '%s\tecoli\t10\tadk(10) fumC(11)\n' "$in"
    ''',
    'mob_recon': '''
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    printf 'sample_id\tmolecule_type\n' > "$out/contig_report.txt"
    printf 'sample_id\tsize\trep_type(s)\trelaxase_type(s)\tpredicted_mobility\n' > "$out/mobtyper_results.txt"
    printf 'asm:plasmid_AA474\t89000\tIncFIB\tMOBF\tconjugative\n' >> "$out/mobtyper_results.txt"
    ''',
    'amrfinder': '''
    [ "$1" = "-u" ] && exit 0
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "-o" ] && out=$a; prev=$a; done
    printf 'Gene symbol\tSequence name\tClass\tSubclass\t%% Coverage of reference sequence\t%% Identity to reference sequence\tContig id\n' > "$out"
    printf 'blaEC\tclass C beta-lactamase\tBETA-LACTAM\tCEPHALOSPORIN\t100.00\t99.75\tcontig_1\n' >> "$out"
    ''',
    'abricate': '''
    db=""; prev=""
    for a in "$@"; do [ "$prev" = "--db" ] && db=$a; prev=$a; done
    printf '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%%COVERAGE\t%%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n'
    printf 'asm\tcontig_1\t10\t500\t+\tgene_%s\t1-490/490\t===\t0/0\t100.00\t98.00\t%s\tACC\tsome product\t\n' "$db" "$db"
    ''',
    'bakta': '''
    out=""; prefix=""; prev=""
    for a in "$@"; do [ "$prev" = "--output" ] && out=$a; [ "$prev" = "--prefix" ] && prefix=$a; prev=$a; done
    mkdir -p "$out"
    { echo "CDSs: 4300"; echo "tRNAs: 86"; echo "rRNAs: 22"; } > "$out/$prefix.txt"
    : > "$out/$prefix.gff3"
    ''',
    'seqkit': '''
    if [ "$1" = "version" ]; then echo "seqkit v2.8.0"; exit 0; fi
    printf 'file\tformat\ttype\tnum_seqs\tsum_len\tN50\tAvgQual\tGC(%%)\n'
    prev=""
    for a in "$@"; do
      case "$a" in stats|-a|-T) continue;; -j) prev=j; continue;; esac
      if [ "$prev" = "j" ]; then prev=""; continue; fi
      printf '%s\tFASTQ\tDNA\t3\t36\t12\t40.0\t50.0\n' "$a"
    done
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
    # A stub database folder so the optional modules that need one (mash species ID) run too
    db_dir = tmp_path_factory.mktemp('stub_db')
    (db_dir / 'mash').mkdir()
    (db_dir / 'mash' / 'refseq.msh').write_bytes(b'stub')
    (db_dir / 'bakta' / 'db-light').mkdir(parents=True)
    (db_dir / 'bakta' / 'db-light' / 'version.json').write_text('{}')
    monkeypatch.setenv('HASHBROWN_DB', str(db_dir))
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
