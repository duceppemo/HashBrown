<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/HashBrown/main/docs/images/hashbrown_logo.svg" alt="HashBrown: long-read, short-read and hybrid bacterial genome assembly" width="560">
</p>

HashBrown assembles bacterial genomes from Nanopore long reads, Illumina (or Ion Torrent) short reads, or
both, and ends every run with QC and a per-sample HTML/PDF report: stage-by-stage read and assembly
funnels, coverage, consensus QV, species confirmation, completeness, typing, AMR genes, annotation and
contamination screens ([Methods](Methods) has the full list; the deeper modules are optional and
auto-detected). Like its namesake — shredded potatoes pressed into one piece — it presses reads into one
assembly. It merges and supersedes [BGA](https://github.com/duceppemo/BGA) and
[BGA_short](https://github.com/duceppemo/BGA_short).

<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/HashBrown/main/docs/images/pipeline.svg" alt="The HashBrown pipeline" width="690">
</p>

## The three ways to run it

```bash
hashbrown -l nanopore_fastq/ -o results/                      # long reads only
hashbrown -l nanopore_fastq/ -s illumina_fastq/ -o results/   # hybrid (short reads polish)
hashbrown -s illumina_fastq/ -o results/                      # short reads only
```

## Pages

- [Installation](Installation) — conda, pip, troubleshooting
- [Usage](Usage) — inputs, options, resuming, performance
- [Methods](Methods) — what each step runs and why, with references
- [Outputs](Outputs) — folders and files
- [Validation](Validation) — benchmarking the polisher combinations against a gold standard
- [FAQ](FAQ) — troubleshooting
- [Development](Development) — tests, CI, releases
