# Changelog

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
