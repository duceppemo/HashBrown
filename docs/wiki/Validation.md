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
