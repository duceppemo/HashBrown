# Reporting & extended QC — design plan (2026-10)

Goal: every run ends with a self-contained HTML report per sample (plus a lightweight PDF twin and a
run-level index comparing samples) that shows **the effect of every pipeline stage** on the reads and on
the assembly, plus species confirmation, completeness, contamination, typing and AMR.

Decisions taken (2026-10-26):

- **Packaging — optional modules.** The always-installed core stays small (seqkit, mosdepth, QUAST,
  pycoQC-ng, jinja2, matplotlib, weasyprint). Heavy tools (CheckM2, BUSCO, Bakta, MOB-suite, merqury,
  sylph, ConFindr, GenomeScope2, AMRFinderPlus, abricate) are **auto-detected at runtime**: installed →
  they run; missing → the report shows "skipped (not installed)" with the install one-liner. Databases are
  fetched by a new `hashbrown db <module>` command. Nothing new may break the main env solve.
- **Read QC engine — pycoQC-ng for both platforms.** Its converters build a sequencing summary from plain
  fastq (`Fastq_to_seq_summary`, recomputing per-read metrics incl. %GC), so it runs on **every stage**:
  raw → trimmed (porechop / fastp) → filtered (filtlong), Nanopore and Illumina alike. A provided
  `sequencing_summary.txt` (new `--summary` flag) unlocks the pore/time plots for raw Nanopore reads.
  FastQC + MultiQC remain as short-read extras; GenomeScope2 (KMC histogram) estimates genome size and
  coverage from short-read k-mers; ConFindr screens for intra-species contamination (full rMLST db needs a
  free PubMLST account — report says so when missing).
- **Annotation — Bakta by default** (when installed; light db), `--annotate pgap` opt-in via the official
  Docker image (~30 GB, ~1 h/genome) for submission-grade output.
- **AMR / typing — AMRFinderPlus + abricate (VFDB, PlasmidFinder)**, outputs normalized into one
  gene × tool table with hAMRonization. MLST via `mlst`; plasmid typing/mobility via MOB-suite.
  (RGI/CARD and ResFinder: not in scope — revisit on demand.)
- **Contamination visual — built-in lightweight blob plot**: GC vs coverage per contig (mosdepth).
  Implementation note (2026-10-26): colored by per-contig `mash dist -i` top hits against the RefSeq
  sketch already used for species ID, named via NCBI datasets — not sylph, which profiles at sample
  level and cannot color individual contigs. No new tool or database needed. Full BlobToolKit stays a
  documented Docker extra.

## Report content

**Read funnel** (one table per read type; one row per stage: raw / trimmed / filtered):
reads, bases, N50, mean Q, %GC, estimated coverage, % of raw surviving — seqkit/nanoq + fastp JSON.
Figures: per-stage length and quality distributions (overlaid), pycoQC-ng per-stage reports linked.

**Assembly funnel** (one row per stage: assembled → medaka → polypolish → pypolca → reoriented):
contigs, size, N50, circular replicons (Flye info), positions changed by each polisher (tool logs),
QUAST metrics, Merqury QV (reference-free consensus quality from short-read k-mers — the in-run version
of the polisher benchmark). Figures: contig coverage plot (mosdepth), dot plots (existing), assembly
graph (existing Bandage image), blob plot.

**Identity & completeness**: mash species ID (GTDB-reps sketch default; proGenomes4 optional via
`hashbrown db`) → auto-fetch the species' type-material assembly (NCBI datasets) → FastANI/skani vs type
strain; CheckM2 completeness/contamination; BUSCO auto-lineage; rRNA/tRNA counts (barrnap + tRNAscan-SE)
as a fast annotation sanity check.

**Typing & genes**: MLST, MOB-suite, AMRFinderPlus + abricate merged table, Bakta/PGAP summary.

**New pipeline step**: dnaapler reorients chromosome (dnaA) and plasmids (rep) after polishing, before
final QC, so every downstream coordinate refers to the fixed orientation.

## Output layout

```
8_report/
  <sample>.html        self-contained (SVG figures inline)
  <sample>.pdf         WeasyPrint render of the same content
  index.html           run-level comparison of all samples
  tables/*.tsv         every table, machine-readable
  pycoqc/, multiqc/    per-tool reports linked from the HTML
```

## Phases

- **P1 — skeleton + funnels (no big DBs):** `hashbrown/report.py` (jinja2 template, matplotlib SVG
  figures, WeasyPrint PDF), read funnel (seqkit + fastp JSON + pycoQC-ng per stage), assembly funnel
  (QUAST multi-assembly + polisher logs), mosdepth per-contig coverage table/plot, run index. Tests with
  stubbed tools, like the existing suite.
- **P2 — identity & completeness:** dnaapler step, mash ID → datasets → FastANI, CheckM2, BUSCO,
  Merqury QV, rRNA/tRNA, `hashbrown db`.
- **P3 — typing, AMR, annotation:** MLST, MOB-suite, AMRFinderPlus, abricate, hAMRonization table,
  Bakta default / PGAP opt-in.
- **P4 — contamination:** ConFindr, sylph, built-in blob plot; GenomeScope2.

Status 2026-10-26: all four phases implemented and pushed (P1 funnels/report, P2 identity &
completeness, P3 typing/AMR/annotation, P4 contamination & GenomeScope2).

Each phase lands as its own commit(s) with tests and wiki updates (Outputs.md gains a "Report" section;
Usage.md the new flags: `--summary`, `--annotate`, `--qc-level`, `hashbrown db`).
