# Outputs

Folders are numbered in pipeline order. A step that was not requested leaves no folder. `log.txt` records
the version of every program used; the `done_*` flag files drive [resuming](Usage#resuming).

## Long reads / hybrid

```
output/
├── log.txt                        program versions
├── flye_stats.tsv                 per sample: read length, N50, assembly length, contigs, coverage
│                                  (shasta_stats.tsv with -a shasta)
├── 1_trimmed_long/                --trim-long (opt-in): Porechop output, <sample>.fastq.gz
├── 2_filtered_long/               Filtlong output (default; absent with --no-filter-long)
├── 3_assembled_long/              <sample>.fasta — the raw assemblies
├── 4_polished_long/               <sample>.fasta — after Medaka
├── 5a_trimmed_short/              fastp output pairs (default; absent with --no-trim-short)
├── 5b_trimmed_short_report/       fastp HTML/JSON reports
├── 6_polished_short/              <sample>.fasta — the FINAL assemblies (contigs renamed contig_N),
│                                  plus <sample>.polypolish.log, <sample>.report and <sample>.vcf (pyPOLCA)
└── 7_assembly_qc/
    ├── quast/                     report.html / report.tsv over all final assemblies
    ├── mapping_qc_long_reads/     <sample>_<stage>.flagstat.txt + <sample>_<stage>.mosdepth.summary.txt
    │                              for EVERY read stage present (raw, trimmed, filtered); the most
    │                              processed stage also under <sample>.flagstat.txt /
    │                              <sample>.mosdepth.summary.txt; long_read_cov.tsv has one row per
    │                              sample and stage
    ├── mapping_qc_short_reads/    the same with the short reads (hybrid runs), short_read_cov.tsv
    ├── mummer_long/               <sample>_medaka.ps dot plot: raw vs Medaka
    ├── mummer_short/              <sample>_polypolish.ps dot plot: Medaka vs short-read polished
    └── assembly_graphs/           <sample>.gfa + <sample>.png (Bandage)
└── 8_report/                     ← examples: a clean run
    │                              ([HTML](https://htmlpreview.github.io/?https://github.com/duceppemo/HashBrown/blob/main/docs/example_report/MG1655.html) ·
    │                              [PDF](https://github.com/duceppemo/HashBrown/blob/main/docs/example_report/MG1655.pdf))
    │                              and a contaminated one
    │                              ([HTML](https://htmlpreview.github.io/?https://github.com/duceppemo/HashBrown/blob/main/docs/example_report/contamdemo.html) ·
    │                              [PDF](https://github.com/duceppemo/HashBrown/blob/main/docs/example_report/contamdemo.pdf))
    ├── index.html                 run overview: one row per sample
    ├── <sample>.html              the sample's report: read funnel (raw → trimmed → filtered),
    │                              assembly stage by stage (with Merqury QV), species confirmation,
    │                              completeness, per-contig coverage, figures, links
    ├── <sample>.pdf               print twin of the HTML (needs WeasyPrint, skipped otherwise)
    ├── tables/                    every report table as TSV
    ├── qc/<sample>/               raw outputs of the optional modules (CheckM2, BUSCO, Merqury, ANI)
    └── pycoqc/                    per-stage pycoQC-ng reports (needs pip install pycoqc-ng)
```

When [dnaapler](https://github.com/gbouras13/dnaapler) is installed, the final assembly of each sample is
reoriented in place (chromosome starts at *dnaA*, plasmids at their *rep* gene) and the pre-reorientation
sequence is kept beside it as `<sample>_unoriented.fasta`; the dot plots always compare pre-reorientation
sequences, so a rotation never shows up as a structural change.

The **final assembly** of each sample is the deepest one present: `6_polished_short/`, else
`4_polished_long/`, else `3_assembled_long/`.

## Short reads only

```
output/
├── log.txt
├── 1_trimmed_short/               fastp output (default; absent with --no-trim-short)
├── 2_assembled_short/             <sample>.fasta — the final assemblies
└── 3_assembly_qc/
    ├── quast/
    ├── mapping_qc_short_reads/    flagstat + mosdepth summaries, short_read_cov.tsv
    └── assembly_graphs/           gfa_connector (SKESA) or SPAdes graphs + Bandage PNGs
└── 4_report/                      same layout as the hybrid 8_report/
```

## File notes

- Every assembly fasta is wrapped at 80 characters; final hybrid assemblies have contigs renamed
  `contig_1`, `contig_2`, … in input order.
- The dot plots are PostScript (`.ps`); any viewer or `ps2pdf` opens them. A straight diagonal means
  polishing changed only small things; steps, gaps or off-diagonal segments mean structural changes worth
  investigating.
- Coverage is mosdepth's overall mean depth (duplicates removed before counting); `NA` means the
  mapping produced no bam (no assembly or no reads).

## What contamination looks like

`docs/example_report/contamdemo` is a deliberately contaminated sample — real *E. coli* K-12 reads
spiked with ~20% *Listeria innocua* reads (both read types), assembled with the defaults. Side by side
with the clean MG1655 report:

| Signal | Clean MG1655 | Contaminated |
| --- | --- | --- |
| MIMAG verdict | High-quality draft | **Low-quality draft** |
| ConFindr | Clean | **Contaminated** ("Escherichia:Listeria") |
| CheckM2 contamination | 0.14% | **66.65%** |
| BUSCO | 99.5% single-copy | **93.1% duplicated** |
| Assembly vs GenomeScope size | 4.64 Mb ≈ 4.66 Mb | **7.61 Mb vs 4.9 Mb** |
| rRNA / tRNA counts | 8/7/7 · 89 | **14/13/13 · 155** (doubled) |
| Merqury QV | ~50 | **~39** |
| Blob plot | one cluster | **two clusters** (E. coli 51%GC/100x, Listeria 37%GC/30x) |

Note what does **not** catch it: the mash species ID still reports *E. coli* at 100% ANI, because the
majority genome dominates — a species check alone is not a contamination check. Several independent
signals firing together is the pattern to trust.
