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
   a 64 bp length floor; adapters auto-detected for pairs. With `--ion-torrent`, adapter trimming is
   disabled (the Torrent Suite removes adapters before export, and fastp's detection is tuned to
   Illumina chemistry); the quality window stays, since Ion quality also decays toward the 3' end.
6. **Polish, short reads**: Polypolish (BWA `-a` all-alignments mapping of each mate, as its author
   specifies), then pyPOLCA with `--careful`. One Polypolish round only: a second round changed nothing in
   our tests. This pair was chosen because neither introduces errors in low-coverage regions —
   NextPolish, ntEdit and plain POLCA did, and were dropped. The final contigs are renamed `contig_1`,
   `contig_2`, … **Polish only with Illumina reads**: Ion Torrent reads share Nanopore's
   homopolymer/indel error mode, so they can confirm those errors instead of fixing them — with
   `--ion-torrent` HashBrown warns and suggests `--short-polisher none`.

## Short reads only

- **Assemble**: SKESA (conservative, recommended for Illumina isolates) or SPAdes (scaffolds;
  `--ion-torrent` adds SPAdes' Ion Torrent mode; single-end supported by both). SKESA graphs are rebuilt
  with `gfa_connector` when it is installed.

## Reorient (every run, optional)

- **dnaapler** (`dnaapler all`), when installed, rotates every final assembly to the conventional start:
  the chromosome at *dnaA*, plasmids at their *rep* gene. The reorientation happens in place and the
  pre-reorientation sequence is kept as `<sample>_unoriented.fasta`. The dot plots always compare
  pre-reorientation sequences, so a rotation never masquerades as a structural change.

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

## Report (every run)

Every run ends with a per-sample HTML report (plus a PDF twin and a run-level index — see
[Outputs](Outputs)). Its core needs nothing beyond the HashBrown environment:

- **Read funnel** per read type: seqkit statistics of every stage (raw → trimmed → filtered), with
  per-stage length/quality distributions and, when `pycoqc-ng` is pip-installed, a full
  [pycoQC-ng](https://github.com/duceppemo/pycoqc-ng) page per stage (its converters rebuild the
  sequencing summary from fastq, for Nanopore and Illumina alike).
- **Assembly funnel**: contigs, size, N50, %GC — and, with Merqury installed, the **consensus QV** of
  every polishing stage from the short-read 21-mers (reference-free; the in-run version of our
  [validation](Validation) measurements).
- **Species confirmation**: mash against a local RefSeq sketch (`hashbrown db mash`), the top hit named
  and fetched with NCBI datasets, ANI with fastANI or skani (same species: ≥95%).
- **Completeness & sanity**: CheckM2 completeness/contamination, BUSCO (prokaryote auto-lineage),
  barrnap rRNA counts, tRNAscan-SE tRNA counts.
- **Typing & genes**: MLST, MOB-suite plasmid reconstruction and typing, and one merged gene table from
  AMRFinderPlus (resistance), abricate+VFDB (virulence) and abricate+PlasmidFinder (replicons).
- **Annotation**: Bakta by default when installed (`--annotate`), NCBI PGAP via Docker on request.
- **Contamination**: ConFindr on the short reads, a GenomeScope2 read-based genome profile, and a
  **blob plot** — every contig by %GC and coverage, dot area ~ length, colored by its closest RefSeq
  genome (per-contig `mash dist -i`).

Each of these modules is optional: a missing tool or database turns into a "skipped (install hint)" note
in the report, never an error. The [Usage](Usage) page lists the install one-liners and `hashbrown db`
downloads. The version of every program used is written to `log.txt` in the output folder.

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
- seqkit — Shen W, et al. SeqKit2: a Swiss army knife for sequence and alignment processing. *iMeta* 3, e191 (2024).
- pycoQC-ng — Duceppe MO, https://github.com/duceppemo/pycoqc-ng (fork of pycoQC: Leger A, Leonardi T. *J Open Source Softw* 4, 1236 (2019)).
- dnaapler — Bouras G, et al. Dnaapler: a tool to reorient circular microbial genomes. *J Open Source Softw* 9, 5968 (2024).
- Merqury — Rhie A, et al. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21, 245 (2020).
- Mash — Ondov BD, et al. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol* 17, 132 (2016).
- FastANI — Jain C, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun* 9, 5114 (2018).
- CheckM2 — Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20, 1203–1212 (2023).
- BUSCO — Manni M, et al. BUSCO update: novel and streamlined workflows. *Mol Biol Evol* 38, 4647–4654 (2021).
- barrnap — Seemann T, https://github.com/tseemann/barrnap
- tRNAscan-SE — Chan PP, et al. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res* 49, 9077–9096 (2021).
- mlst — Seemann T, https://github.com/tseemann/mlst (uses PubMLST: Jolley KA, et al. *Wellcome Open Res* 3, 124 (2018)).
- MOB-suite — Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb Genom* 4, e000206 (2018).
- AMRFinderPlus — Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog. *Sci Rep* 11, 12728 (2021).
- abricate — Seemann T, https://github.com/tseemann/abricate (VFDB: Liu B, et al. *Nucleic Acids Res* 50, D912–D917 (2022); PlasmidFinder: Carattoli A, et al. *Antimicrob Agents Chemother* 58, 3895–3903 (2014)).
- Bakta — Schwengers O, et al. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 7, 000685 (2021).
- PGAP — Tatusova T, et al. NCBI prokaryotic genome annotation pipeline. *Nucleic Acids Res* 44, 6614–6624 (2016).
- ConFindr — Low AJ, et al. ConFindr: rapid detection of intraspecies and cross-species contamination in bacterial whole-genome sequence data. *PeerJ* 7, e6995 (2019).
- GenomeScope2 — Ranallo-Benavidez TR, et al. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nat Commun* 11, 1432 (2020).
- WeasyPrint — Kozea, https://weasyprint.org
