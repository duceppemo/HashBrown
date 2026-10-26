# Changelog

## Unreleased

### Changed
- **Single-end short reads can polish long-read assemblies**: Polypolish (one all-alignments SAM) and
  pyPOLCA (no `--reads2`) both run single-end, so Ion Torrent or unpaired Illumina data works in hybrid
  mode; a sample with a pair plus stray extra files keeps the pair with a warning. New `--skip-report`
  flag for harness runs where only the assemblies matter.
- **`--ion-torrent` now shapes the whole run, not just SPAdes**: fastp skips adapter trimming for Ion
  Torrent reads (the Torrent Suite already removed the adapters; the sliding-window quality trim stays),
  the flag is accepted in hybrid mode, and polishing a Nanopore assembly with Ion Torrent reads prints a
  warning — both platforms share homopolymer/indel errors, so such a polish can confirm errors rather
  than fix them (`--short-polisher none` is suggested).

### Added
- **Contamination screening and read-based profiling** (reporting phase 4): ConFindr verdict from the
  short reads (`hashbrown db confindr`), a GenomeScope2 genome profile (size estimate, heterozygosity,
  read error rate from the short-read 21-mers, reusing Merqury's meryl counts), and a built-in **blob
  plot** — every contig by %GC and coverage, dot area ~ length, colored by its closest RefSeq genome
  (per-contig `mash dist -i` + NCBI datasets; no BlobToolKit-scale databases needed). The plot plus a
  `tables/<sample>_blob.tsv` make stray contigs obvious.
- **Typing, AMR and annotation in the report** (reporting phase 3), optional modules like the rest:
  MLST (scheme/ST also on the run index), MOB-suite plasmid reconstruction and typing, one merged
  gene table from AMRFinderPlus (resistance), abricate+VFDB (virulence) and abricate+PlasmidFinder
  (replicons), and genome annotation — Bakta by default when installed, NCBI PGAP via Docker with
  `--annotate pgap` (new flag: `--annotate {auto,bakta,pgap,none}`). `hashbrown db` gains the
  `bakta`, `amrfinder` and `mob` modules.
- **Optional QC modules** (reporting phase 2), auto-detected at runtime and reported as skipped-with-a-hint
  when not installed: **dnaapler** reorients every final assembly in place (chromosome to *dnaA*, plasmids
  to *rep*; the input order is kept next to it as `<sample>_unoriented.fasta`, and dot plots keep comparing
  pre-reorientation sequences); **Merqury** consensus QV per assembly stage from the short-read k-mers (new
  QV column in the report's assembly table); **species confirmation** — mash against a local RefSeq sketch,
  then ANI (fastANI or skani) against the top hit fetched with NCBI `datasets`; **CheckM2**
  completeness/contamination; **BUSCO** (prokaryote auto-lineage); **barrnap** rRNA counts and
  **tRNAscan-SE** tRNA counts as annotation sanity checks. Results land in the HTML/PDF report, the
  `tables/` TSVs and `8_report/qc/<sample>/`.
- **`hashbrown db {mash,checkm2}`** downloads the databases those modules need (default
  `~/.hashbrown/db`, override with `--dir` or `$HASHBROWN_DB`).
- **Per-sample HTML report** (`8_report/`, or `4_report/` for short-read-only runs) with a run-level
  `index.html`: read funnel tables showing the effect of every read-processing stage (raw → trimmed →
  filtered: reads, bases, N50, mean Q, %GC, coverage, % kept), assembly stage-by-stage table
  (assembled → Medaka → short-read polished), per-contig mosdepth coverage (long and short reads side by
  side), length/quality distribution and coverage figures, and links to the QUAST/pycoQC/dot plot details.
  Every table is also written as TSV under `8_report/tables/`. A PDF twin is written when WeasyPrint is
  installed, and per-stage [pycoQC-ng](https://github.com/duceppemo/pycoqc-ng) reports when
  `pip install pycoqc-ng` was run — both are optional and skipped with a note otherwise.
  New core dependencies: seqkit, jinja2, matplotlib. The full reporting/QC roadmap is in
  `docs/reporting_plan.md`.

## 1.0.0 (2026-10-24)

HashBrown is the merge of [BGA](https://github.com/duceppemo/BGA) (long-read and hybrid assembly) and
[BGA_short](https://github.com/duceppemo/BGA_short) (short-read assembly) into one tool with one conda
environment. The split existed because Medaka's TensorFlow builds conflicted with the short-read tools;
Medaka ≥ 2.0 is PyTorch-based and the conflict is gone. Both repositories' histories are kept in this one.

### Changed
- **One pipeline for every input**: `-l` long reads only, `-l` + `-s` hybrid (the short reads polish the
  long-read assemblies), or `-s` short reads only (SKESA by default, or SPAdes; single-end or paired-end;
  `--ion-torrent` for SPAdes).
- Renamed and repackaged: `pip install .` gives the `hashbrown` command (package `hashbrown`, Python ≥ 3.10);
  `python hashbrown.py` still works from a clone. Pandas is no longer needed; `psutil` is the only Python
  dependency.
- **The polishers are now chosen on the command line**, so combinations can be benchmarked against a gold
  standard (`validation/benchmark_polishers.sh`): `--long-polisher {medaka,none}` and
  `--short-polisher {polypolish+pypolca,polypolish,pypolca,none}`. The defaults are unchanged: Medaka, then
  Polypolish followed by pyPOLCA `--careful` when short reads are given — the combination that introduced no
  errors in low-coverage regions in earlier tests. Short-read polishing now happens whenever short reads
  are provided alongside long reads (the old `--polish` flag is gone); pass `--short-polisher none` to
  turn it off.
- `--show-models` is now `--list-models`, and no longer hangs when medaka's help does not list the models.
- Dropped tools that were no longer called: NextPolish, ntEdit, BBMap, Falco, pycoQC and the second and
  third Polypolish rounds (one round changes nothing further; some of the removed polishers introduced
  errors in low-coverage regions). fastp's reports cover the short-read QC.
- The SKESA assembly graphs use `gfa_connector` from the bioconda SKESA package when it is on the PATH
  instead of a bundled binary; Bandage pictures and gfa_connector are optional (a missing program skips
  that output with a note instead of failing).
- **Qualimap is replaced by `samtools flagstat` + mosdepth**: same information (mapping QC and depth),
  one mapping pass instead of two, and it removes the Java/R stack that made the conda environment
  unsolvable. The QC folders are `mapping_qc_long_reads/` and `mapping_qc_short_reads/`, holding the
  flagstat and mosdepth summaries and the coverage tables.

- **The default Medaka model is `r1041_e82_400bps_bacterial_methylation`** (ONT's bacterial model):
  in the gold-standard benchmarks (see the wiki's Validation page) it removed 3x more indels than the
  matched sup model on real R10.4.1 data and never performed worse. The defaults
  `--long-polisher medaka --short-polisher polypolish+pypolca` are confirmed by the same benchmarks;
  use `--long-polisher none` for non-Nanopore long reads (Medaka on PacBio doubles the mismatch rate).

### Fixed
- The thread/parallel/memory checks actually clamp the values now (the old comparisons could never be true,
  so invalid values passed through).
- Hybrid mode mapped the short reads into the long-read Qualimap folder, so the short-read Qualimap reports
  were computed from nothing.
- The polished-vs-polished dot plots are named `<sample>_polypolish.ps`; they used to overwrite the
  `<sample>_medaka.ps` plots of the raw-vs-medaka comparison.
- The coverage tables are written once per run instead of appended, so resuming no longer duplicates rows,
  and a sample with an empty bam reports 0 instead of crashing on a division by zero.
- Resuming after the assembly step no longer overwrites `flye_stats.tsv`/`shasta_stats.tsv` with an empty
  table.
- Short-read samples without matching long reads, files where R1/R2 cannot be told apart, and duplicate
  long-read files for one sample are reported and skipped instead of crashing (`KeyError`) or being
  silently mis-paired.
- SPAdes single-end assemblies pass `-s <file>` (the old `--s 1` was not a valid option, so single-end
  SPAdes runs always failed); `--ion-torrent` is actually forwarded to SPAdes now.
- fastp writes its JSON report next to the HTML one instead of leaving `fastp.json` in the working
  directory, which crashed parallel runs that shared it.
- A missing reference or read file no longer leaves minimap2/samtools zombie pipes behind.

---

The history of BGA and BGA_short before the merge is in `git log` (both repositories' commits are kept).
