# Validation: choosing the polishers with data

Polishing is where assemblies silently get *worse*: several short-read polishers used in BGA's past
(NextPolish, ntEdit, POLCA without `--careful`, extra Polypolish rounds) introduced errors in low-coverage
regions. HashBrown therefore exposes the polishers on the command line and ships a benchmark harness so the
defaults follow measurements, not habit:

- `--long-polisher {medaka,none}` — with modern super-accuracy basecalls, Medaka's benefit shrinks; measure
  it on your own data before trusting either choice.
- `--short-polisher {polypolish+pypolca,polypolish,pypolca,none}` — the default pair is the combination
  that introduced no errors in earlier tests (Polypolish and pyPOLCA `--careful` are both designed to
  leave low-confidence regions alone).

## Benchmark results (2026-10, HashBrown 1.0.0)

Five tracks, every assembly scored by QUAST against a finished genome (mismatches and indels per 100 kbp,
lower is better). Full tables and runs: `benchmark_results.tsv` under each run folder.

**Simulated M. bovis AF2122/97** (badread Q20 Nanopore 100x + ART Illumina 60x; exact truth): Flye alone
left one indel in 4.35 Mb; every polisher combination removed it and none introduced an error. All
combinations are safe at good coverage.

**Real E. coli K-12 MG1655** (MicroBench R10.4.1 PromethION 100x + NovaSeq 60x; truth U00096.3): the Flye
assembly was already at Q63 (2 mismatch + 6 indel positions, which are real lab-stock differences from
U00096.3). Medaka (default and bacterial model) and the short-read polishers changed at most 2 positions
and nothing got worse. With recent super-accuracy basecalls at 100x, polishing has almost nothing left
to fix.

**Real M. bovis field strain** (R10.4.1 sup v4.3.0, 100x, compared against AF2122/97 — strain divergence
is a constant offset, so only the differences between arms count):

| Long polisher | Mismatches/100kbp | Indels/100kbp |
| --- | --- | --- |
| none | 19.12 | 5.42 |
| medaka, matched sup model | 19.10 | 5.26 |
| medaka, bacterial model | 19.10 | **4.64** |

Medaka never hurt, and **ONT's bacterial methylation model removed 3x more indels than the matched sup
model** (~34 vs ~7 genome-wide). It is HashBrown's default model since this benchmark.

**Real AF2122/97 PacBio RS CLR + NextSeq** (truth = its own finished genome; the CLR base assembly is
deliberately rough — a stress test):

| Combination | Mismatches/100kbp | Indels/100kbp |
| --- | --- | --- |
| none + none | 52 | 516 |
| none + polypolish | 45 | 362 |
| none + pypolca | 52 | 511 |
| none + polypolish→pypolca | **45** | **361** |
| medaka + none | 100 | 594 |

Polypolish does the heavy lifting and never hurts; pyPOLCA `--careful` alone declines to fix a bad
assembly (by design) but adds a small gain after Polypolish. **Medaka on PacBio reads doubles the
mismatches: never run Medaka on non-Nanopore reads** (`--long-polisher none`).

**Real MG1655 PacBio HiFi + NovaSeq** (Flye `--pacbio-hifi` base): the HiFi assembly was near-perfect and
short-read polishing made it marginally worse (+1 mismatch, +2 indels) because the Illumina reads come
from a different lab stock. Polishing an already-excellent assembly with reads from a different culture
moves it toward that culture.

### What the defaults follow from this

- **Hybrid Nanopore + Illumina**: keep `--long-polisher medaka --short-polisher polypolish+pypolca`
  (the default). Gains are small on good data, large on rough assemblies, and never negative on Nanopore.
- **Medaka model**: `r1041_e82_400bps_bacterial_methylation` (the default) for bacterial isolates on
  R10.4.1; it beat or equalled the matched sup model everywhere.
- **Recent sup basecalls at ≥100x**: `--long-polisher none` is acceptable — the assembler is already at
  the error floor — but medaka costs little and still helps older basecalls (sup v4.3: −34 indels).
- **Non-Nanopore long reads** (PacBio): always `--long-polisher none`.
- Polish with reads **from the same culture** as the assembly; reads from a different stock of the same
  strain introduce that stock's variants.

## Running the benchmark

```bash
conda activate HashBrown
bash validation/benchmark_polishers.sh data/ results/ 16
```

`data/` holds one folder per sample: `reference.fasta` (the gold standard), the Nanopore fastq and the
Illumina pair — see [`validation/README.md`](https://github.com/duceppemo/HashBrown/blob/master/validation/README.md)
for the exact layout and how to simulate reads (badread + wgsim) when you have no finished genome of your
own.

The harness assembles each sample **once**, then resumes a copy of that state for each of the 8 polisher
combinations, compares every final assembly with the gold standard using QUAST, and writes
`benchmark_results.tsv`: mismatches and indels per 100 kbp, misassemblies, NGA50, length and contig count
per combination.

## Reading the results

- The best combination has the fewest mismatches + indels per 100 kbp and no misassemblies.
- Indels are the Nanopore signature error; if `none` → `polypolish+pypolca` barely beats
  `medaka` → `polypolish+pypolca`, Medaka is not paying for its runtime on your basecalls.
- Check the low-coverage samples separately: that is where polishers turn harmful, and a mean over
  high-coverage samples hides it.

If your measurements disagree with HashBrown's defaults, open an issue with the table — the defaults follow
the evidence.
