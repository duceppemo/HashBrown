# Methods

What each step runs, with the exact settings that are not defaults.

## Long reads

1. **Trim** (`--trim-long`): Porechop removes adapters (`--check_reads 1000`).
2. **Filter** (`--filter-long`): Filtlong keeps the best 95% of the bases (`--keep_percent 95`,
   `--min_length 500`); with `--size`, at most the best 100x of the genome (`--target_bases`).
3. **Assemble**: Flye (`--nano-hq`, 3 polishing iterations; expects recent super-accuracy basecalls) or
   Shasta (`--config Nanopore-May2022`, minimum read length 3,000 unless `--min-size`). The assembly
   graph is kept and drawn with Bandage; per-sample stats land in `flye_stats.tsv` / `shasta_stats.tsv`.
4. **Polish, long reads**: Medaka (`medaka_consensus`) with the model from `--model`. The default is
   ONT's bacterial methylation model, the best performer in our [benchmarks](Validation); for
   non-bacterial samples or other flowcells pick the model matching your basecaller — a wrong model
   *adds* errors. `--long-polisher none` skips this: do that for recent super-accuracy basecalls at high
   coverage (little left to fix), and **always for non-Nanopore long reads** (Medaka on PacBio doubled
   the mismatch rate in our benchmark).

## Short reads (hybrid polishing)

5. **Trim** (`--trim-short`): fastp with a sliding-window quality cut (`--cut_right`, mean quality 10) and
   a 64 bp length floor; adapters auto-detected for pairs.
6. **Polish, short reads**: Polypolish (BWA `-a` all-alignments mapping of each mate, as its author
   specifies), then pyPOLCA with `--careful`. One Polypolish round only: a second round changed nothing in
   our tests. This pair was chosen because neither introduces errors in low-coverage regions —
   NextPolish, ntEdit and plain POLCA did, and were dropped. The final contigs are renamed `contig_1`,
   `contig_2`, …

## Short reads only

- **Assemble**: SKESA (conservative, recommended for Illumina isolates) or SPAdes (scaffolds;
  `--ion-torrent` adds SPAdes' Ion Torrent mode; single-end supported by both). SKESA graphs are rebuilt
  with `gfa_connector` when it is installed.

## QC (every run)

- **QUAST** over the final assemblies (`--min-contig 1000`).
- **Mapping QC** per sample and read type, on a minimap2 mapping (`map-ont` or `sr` preset, unmapped
  reads dropped, duplicates removed with `samtools markdup`): `samtools flagstat` for the mapping rates
  and mosdepth for the depth per contig.
- **Coverage**: mosdepth's mean depth over the final assembly, one table per read type.
- **Dot plots**: nucmer + mummerplot comparing the assembly before and after each polishing stage
  (raw vs Medaka, Medaka vs short-read polished) — a quick visual check that polishing did not
  rearrange anything.
- **Assembly graphs**: Bandage PNGs of the Flye/Shasta/SPAdes/gfa_connector graphs.

The version of every program used is written to `log.txt` in the output folder.

## References

- Flye — Kolmogorov M, et al. Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol* 37, 540–546 (2019).
- Shasta — Shafin K, et al. Nanopore sequencing and the Shasta toolkit enable efficient de novo assembly of eleven human genomes. *Nat Biotechnol* 38, 1044–1053 (2020).
- Medaka — Oxford Nanopore Technologies, https://github.com/nanoporetech/medaka
- Porechop — Wick RR, https://github.com/rrwick/Porechop
- Filtlong — Wick RR, https://github.com/rrwick/Filtlong
- SKESA (and gfa_connector) — Souvorov A, et al. SKESA: strategic k-mer extension for scrupulous assemblies. *Genome Biol* 19, 153 (2018). https://github.com/ncbi/SKESA
- SPAdes — Prjibelski A, et al. Using SPAdes de novo assembler. *Curr Protoc Bioinformatics* 70, e102 (2020).
- fastp — Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34, i884–i890 (2018).
- Polypolish — Wick RR, Holt KE. Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Comput Biol* 18, e1009802 (2022).
- pyPOLCA — Bouras G, et al. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. *Microb Genom* 10, 001254 (2024). POLCA: Zimin AV, Salzberg SL. *PLoS Comput Biol* 16, e1007981 (2020).
- BWA — Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics* 25, 1754–1760 (2009).
- minimap2 — Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34, 3094–3100 (2018).
- SAMtools — Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10, giab008 (2021).
- QUAST — Gurevich A, et al. QUAST: quality assessment tool for genome assemblies. *Bioinformatics* 29, 1072–1075 (2013).
- mosdepth — Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. *Bioinformatics* 34, 867–868 (2018).
- MUMmer4 — Marçais G, et al. MUMmer4: a fast and versatile genome alignment system. *PLoS Comput Biol* 14, e1005944 (2018).
- Bandage — Wick RR, et al. Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics* 31, 3350–3352 (2015).
