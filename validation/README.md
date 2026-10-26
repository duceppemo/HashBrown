# Validation: which polisher combination is best?

Some short-read polishers used in the past (NextPolish, ntEdit, POLCA without `--careful`, extra Polypolish
rounds) **introduced** errors in low-coverage regions instead of fixing them. HashBrown therefore exposes the
polishers on the command line (`--long-polisher`, `--short-polisher`) and ships this harness to measure every
combination against a gold standard, so the defaults are chosen by data, not habit.

The exact scripts, accessions and gold standards behind HashBrown's published results (the wiki's
Validation page) live in [`benchmarks/`](benchmarks/): six reproducible tracks over two *M. bovis* strains,
*E. coli* K-12 and two non-model isolates, on Nanopore, PacBio and Illumina.

## Running the benchmark

```bash
conda activate HashBrown
bash validation/benchmark_polishers.sh <data_folder> <results_folder> [threads]
```

`<data_folder>` holds one subfolder per sample with `reference.fasta` (the gold standard), the Nanopore reads
and the Illumina pair (see the header of the script for the exact names). The script assembles each sample
once, resumes a copy of that state for each of the 8 combinations (`none`/`medaka` ×
`none`/`polypolish`/`pypolca`/`polypolish+pypolca`), compares every final assembly with the reference using
QUAST, and writes `benchmark_results.tsv`:

| Column | Meaning |
| --- | --- |
| MismatchesPer100kbp | substitution errors against the gold standard (lower is better) |
| IndelsPer100kbp | insertion/deletion errors — the dominant Nanopore error mode (lower is better) |
| Misassemblies | structural differences (should be 0 on a good assembly) |
| NGA50 | contiguity once broken at misassemblies (higher is better) |

## Choosing the gold standard

- **Real data**: a strain with a trusted finished genome (e.g. a completed reference sequenced again in your
  lab). Mismatches then include true biological differences, so compare combinations against each other
  rather than expecting zero.
- **Simulated data**: simulate reads from any finished genome and the reference is the exact truth, e.g.
  `badread simulate --reference ref.fasta --quantity 100x` for Nanopore and
  `wgsim -1 150 -2 150 -e 0.001 ref.fasta R1.fq R2.fq` for Illumina (both on bioconda; not part of the
  HashBrown environment).

## Interpreting

The combination with the fewest mismatches + indels per 100 kbp and no misassemblies wins. Watch the
low-coverage regions: a polisher that helps at 100x can still hurt at 20x, which is exactly how the dropped
polishers failed. Rerun the benchmark whenever a polisher, an assembler or a basecaller version changes, and
report results in an issue or pull request so the defaults can follow the evidence.
