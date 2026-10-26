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
| `--model` | `r1041_e82_400bps_bacterial_methylation` | Medaka model. The default is ONT's bacterial model (R10.4.1, basecaller-independent), the [benchmark](Validation) winner; `--list-models` prints the alternatives |
| `--short-polisher` | `polypolish+pypolca` | `polypolish+pypolca`, `polypolish`, `pypolca` or `none`, applied in order |
| `--ion-torrent` | off | the short reads are Ion Torrent data: fastp skips adapter trimming, SPAdes uses its Ion Torrent mode. Polishing a Nanopore assembly with them is warned against (shared homopolymer errors) |
| `--annotate` | `auto` | report annotation: `auto` (Bakta when installed), `bakta`, `pgap` (Docker, ~1 h/genome) or `none` |
| `-t, --threads` | all | total threads |
| `-p, --parallel` | 2 | samples processed at once (each gets `threads / parallel`) |
| `-m, --memory` | 85% of RAM | memory in GB (SKESA, SPAdes) |

Pick the Medaka model that matches your flowcell and basecaller: `hashbrown --list-models`.

## Optional QC modules

Every run ends with an HTML report (`8_report/`, or `4_report/` without long reads). The core of it —
stage funnels, coverage, figures — needs nothing extra. The deeper QC runs **only when its tool is
installed** (any missing one is reported as "skipped" with the install hint, never an error):

| Tool | Adds to the report | Install |
| --- | --- | --- |
| [pycoQC-ng](https://github.com/duceppemo/pycoqc-ng) | per-stage read QC pages, both platforms | `pip install pycoqc-ng` |
| WeasyPrint | the PDF twin of the report | `conda install weasyprint` |
| dnaapler | final assemblies reoriented to *dnaA*/*rep* | `conda install dnaapler` |
| Merqury | consensus QV per assembly stage (short-read k-mers) | `conda install merqury` |
| mash + NCBI datasets + fastANI (or skani) | species confirmation: closest genome and ANI | `conda install mash ncbi-datasets-cli fastani` |
| CheckM2 | completeness / contamination | `conda install checkm2` |
| BUSCO | single-copy ortholog completeness (auto-lineage) | `conda install busco` |
| barrnap, tRNAscan-SE | rRNA / tRNA sanity counts | `conda install barrnap trnascan-se` |
| mlst | 7-gene MLST scheme and sequence type | `conda install mlst` |
| MOB-suite | plasmid reconstruction, replicon/relaxase typing, mobility | `conda install mob_suite` |
| AMRFinderPlus | resistance genes and point mutations | `conda install ncbi-amrfinderplus` |
| abricate | virulence factors (VFDB) and plasmid replicons (PlasmidFinder) | `conda install abricate` |
| Bakta | genome annotation (default when installed; `--annotate none` to skip) | `conda install bakta` |
| PGAP | submission-grade NCBI annotation via Docker (`--annotate pgap`, ~1 h/genome) | [ncbi/pgap](https://github.com/ncbi/pgap) |
| ConFindr | intra-species contamination from the short reads | `conda install confindr` |
| GenomeScope2 | read-based genome size/heterozygosity estimate (needs merqury's meryl) | `conda install genomescope2` |

The **blob plot** (per-contig GC vs coverage, colored by the closest RefSeq genome) needs no extra tool:
it uses mash, seqkit and NCBI datasets, and draws uncolored when the mash database or internet is absent.

Some modules need a local database first (default `~/.hashbrown/db`, override with `$HASHBROWN_DB`;
AMRFinderPlus and MOB-suite manage their own locations):

```bash
hashbrown db mash       # RefSeq genome sketches for the species ID (~800 MB)
hashbrown db checkm2    # CheckM2 DIAMOND database (~3 GB)
hashbrown db bakta      # Bakta light database (~2 GB)
hashbrown db amrfinder  # AMRFinderPlus database
hashbrown db mob        # MOB-suite database
hashbrown db confindr   # ConFindr rMLST core databases (full coverage needs a free PubMLST account)
```

Install the heavy ones in their own conda env if they ever conflict with HashBrown's, or just to keep the
solve fast — the pipeline only needs them on the `PATH`. CheckM2 and BUSCO pull large dependency stacks,
so this is the recommended way to add them:

```bash
conda create -n hashbrown-qc -c conda-forge -c bioconda checkm2 busco
# expose them to HashBrown runs (conda entry points carry their own env in the shebang):
ln -s "$(conda info --base)/envs/hashbrown-qc/bin/"{checkm2,busco} "$CONDA_PREFIX/bin/"
```

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
