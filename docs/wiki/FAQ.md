# FAQ

**My sample is missing from the results.**
The sample name is everything before the first underscore of the file name. Two files whose names start the
same way collapse into one sample (the extra long-read file is reported and ignored); a long-read file
named `unclassified*` is dropped on purpose (unbinned barcode reads). In hybrid mode, short-read files whose
name matches no long-read sample are reported and skipped.

**Why did short-read polishing not run?**
It needs a *pair* (`_R1` + `_R2`) and an assembly for that sample. Single-end short reads polish nothing
(Polypolish and pyPOLCA want pairs); they are only assembled in short-read-only mode. Also check
`--short-polisher` was not `none`.

**Can I assemble PacBio reads?**
The long-read pipeline is tuned for Nanopore (Flye `--nano-hq`, Medaka). PacBio HiFi will assemble, but use
`--long-polisher none` (Medaka models are Nanopore-only) — or better, a HiFi-specific assembler.

**Medaka crashes or hangs.**
Medaka runs one sample at a time by design (parallel instances crash). If a single run still fails, the
usual culprit is a broken samtools/bcftools in the environment: run `samtools version`, and recreate the
environment if it errors. Also check the `--model` matches your basecaller — `hashbrown --list-models`.

**How do I rerun just one step?**
Delete that step's `done_*` flag and its output folder, then rerun the same command — see
[Usage § Resuming](Usage#resuming).

**It says a program is not installed, but I don't use that step.**
HashBrown only checks the programs your flags need. If it names one, some requested step uses it —
e.g. `nucmer`/`mummerplot`/`gnuplot` belong to the long-read QC, `bwa`/`polypolish`/`pypolca` to hybrid
polishing.

**Out of memory.**
Lower `-p` (samples in parallel) before lowering `-m`: each parallel sample multiplies the footprint.
SPAdes splits `-m` across the parallel samples; SKESA takes it per sample.

**Where is the final assembly?**
The deepest folder that exists for the sample: `6_polished_short/`, else `4_polished_long/`, else
`3_assembled_long/` (short-read-only runs: `2_assembled_short/`). See [Outputs](Outputs).

**How is this different from BGA / BGA_short?**
One tool, one environment, pluggable polishers, plus the fixes listed in the
[changelog](https://github.com/duceppemo/HashBrown/blob/master/CHANGELOG.md). The old repositories stay up
but are superseded.
