# Outputs

Folders are numbered in pipeline order. A step that was not requested leaves no folder. `log.txt` records
the version of every program used; the `done_*` flag files drive [resuming](Usage#resuming).

## Long reads / hybrid

```
output/
├── log.txt                        program versions
├── flye_stats.tsv                 per sample: read length, N50, assembly length, contigs, coverage
│                                  (shasta_stats.tsv with -a shasta)
├── 1_trimmed_long/                --trim-long: Porechop output, <sample>.fastq.gz
├── 2_filtered_long/               --filter-long: Filtlong output
├── 3_assembled_long/              <sample>.fasta — the raw assemblies
├── 4_polished_long/               <sample>.fasta — after Medaka
├── 5a_trimmed_short/              --trim-short: fastp output pairs
├── 5b_trimmed_short_report/       fastp HTML/JSON reports
├── 6_polished_short/              <sample>.fasta — the FINAL assemblies (contigs renamed contig_N),
│                                  plus <sample>.polypolish.log, <sample>.report and <sample>.vcf (pyPOLCA)
└── 7_assembly_qc/
    ├── quast/                     report.html / report.tsv over all final assemblies
    ├── mapping_qc_long_reads/     <sample>.flagstat.txt + <sample>.mosdepth.summary.txt
    │                              and long_read_cov.tsv (mean depth per sample)
    ├── mapping_qc_short_reads/    the same with the short reads (hybrid runs), short_read_cov.tsv
    ├── mummer_long/               <sample>_medaka.ps dot plot: raw vs Medaka
    ├── mummer_short/              <sample>_polypolish.ps dot plot: Medaka vs short-read polished
    └── assembly_graphs/           <sample>.gfa + <sample>.png (Bandage)
```

The **final assembly** of each sample is the deepest one present: `6_polished_short/`, else
`4_polished_long/`, else `3_assembled_long/`.

## Short reads only

```
output/
├── log.txt
├── 1_trimmed_short/               --trim-short: fastp output (report/ holds the HTML/JSON)
├── 2_assembled_short/             <sample>.fasta — the final assemblies
└── 3_assembly_qc/
    ├── quast/
    ├── mapping_qc_short_reads/    flagstat + mosdepth summaries, short_read_cov.tsv
    └── assembly_graphs/           gfa_connector (SKESA) or SPAdes graphs + Bandage PNGs
```

## File notes

- Every assembly fasta is wrapped at 80 characters; final hybrid assemblies have contigs renamed
  `contig_1`, `contig_2`, … in input order.
- The dot plots are PostScript (`.ps`); any viewer or `ps2pdf` opens them. A straight diagonal means
  polishing changed only small things; steps, gaps or off-diagonal segments mean structural changes worth
  investigating.
- Coverage is mosdepth's overall mean depth (duplicates removed before counting); `NA` means the
  mapping produced no bam (no assembly or no reads).
