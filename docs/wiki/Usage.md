# Usage

```
hashbrown [-l LONG_FOLDER] [-s SHORT_FOLDER] -o OUTPUT_FOLDER [options]
```

Give long reads (`-l`), short reads (`-s`), or both; the pipeline follows the inputs:

| Inputs | Pipeline |
| --- | --- |
| `-l` | long-read assembly (Flye or Shasta) + Medaka polishing + QC |
| `-l` + `-s` | the same, then short-read polishing (Polypolish → pyPOLCA) |
| `-s` | short-read assembly (SKESA or SPAdes) + QC |

## Inputs

- **Long reads**: one fastq file per sample (gzipped or not, found recursively). The sample name is
  everything before the first underscore of the file name — `ecoli_pass.fastq.gz` is sample `ecoli`.
  A file named `unclassified*` (Nanopore barcoding leftovers) is dropped.
- **Short reads**: `_R1`/`_R2` in the file name mark a pair; a file without either is a single-end sample
  (short-read-only pipeline; the polishers need pairs, so in hybrid mode such files are skipped).
- **Hybrid**: the short-read sample names must match the long-read ones (same text before the first
  underscore). Unmatched short-read files are reported and skipped.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `-a, --assembler` | `flye` with `-l`, `skesa` with `-s` alone | `flye`, `shasta` (long reads), `skesa`, `spades` (short reads only) |
| `--trim-long` | off | trim long-read adapters with Porechop before assembling |
| `--filter-long` | off | drop the worst 5% of the long reads (and reads < 500 bp) with Filtlong |
| `--trim-short` | off | trim the short reads with fastp (quality ≥ 10 sliding window, length ≥ 64) |
| `--size` | auto | genome size for Flye; also caps the filtered reads at the best 100x |
| `--min-size` | Shasta 3000 / Flye auto | minimum read length (Shasta) or minimum overlap (Flye) |
| `--long-polisher` | `medaka` | `medaka` or `none` |
| `--model` | `r1041_e82_400bps_sup_g615` | Medaka model; `--list-models` prints the available ones |
| `--short-polisher` | `polypolish+pypolca` | `polypolish+pypolca`, `polypolish`, `pypolca` or `none`, applied in order |
| `--ion-torrent` | off | the short reads are Ion Torrent data (SPAdes only) |
| `-t, --threads` | all | total threads |
| `-p, --parallel` | 2 | samples processed at once (each gets `threads / parallel`) |
| `-m, --memory` | 85% of RAM | memory in GB (SKESA, SPAdes) |

Pick the Medaka model that matches your flowcell and basecaller: `hashbrown --list-models`.

## Resuming

Each major step leaves a `done_*` flag file in the output folder. Rerunning the same command skips the
steps already done and continues; the QC always reruns (it is cheap and depends on everything before it).
To redo a step, delete its flag (and its output folder, e.g. `done_polishing_short` and
`6_polished_short/`) and rerun.

## Performance

- `-p/--parallel` trades memory for speed: each of the `p` samples runs with `threads / p` threads. Flye
  and SPAdes can use tens of GB per sample — keep `-p 2` (or 1) unless the machine is large.
- Medaka always runs one sample at a time with all the threads: parallel Medaka instances crash.
- The QC maps reads once per report; with many samples it is the longest part after assembly.
