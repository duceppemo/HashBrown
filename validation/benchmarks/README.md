# Benchmark tracks (2026-10)

These are the exact scripts behind the results in the wiki's
[Validation](https://github.com/duceppemo/HashBrown/wiki/Validation) page. Each track feeds
[`../benchmark_polishers.sh`](../benchmark_polishers.sh) (or a small arm-loop of its own) and scores every
final assembly against a gold-standard genome with QUAST (`-r`). Activate the HashBrown conda env before
running any of them; `badread`, `art` and `seqtk` are only needed for data preparation and are not part of
the HashBrown environment.

| Track | Script | Organism | Reads | Truth | Public data |
| --- | --- | --- | --- | --- | --- |
| AFsim | `afsim_make_data.sh` + harness | *M. bovis* AF2122/97 | simulated ONT (badread Q20) + ART Illumina | exact (simulation source LT708304.1) | yes |
| AF2122 | `af2122_benchmark.sh` | *M. bovis* AF2122/97 | real PacBio RS CLR (ERR1748639) + NextSeq (ERR4586795) | its own finished genome LT708304.1 | yes |
| MG1655 | `mg1655_benchmark.sh` | *E. coli* K-12 (DSM 18039) | real ONT R10.4.1 (ERR17025883) + Illumina (SRR41214148) + HiFi (SRR10971019) | U00096.3 (lab-stock caveat) | yes |
| MBWGS426 | `mbwgs426_medaka_models.sh` | *M. bovis* field strain | real ONT R10.4.1 sup v4.3.0 (in-house) | AF2122/97 (comparative deltas only) | reads in-house |
| Non-model | `autocycler_nonmodel_benchmark.sh` | *L. innocua*, *P. rettgeri* | real ONT R10.4.1 sup v5.0.0 + Illumina | same-isolate Trycycler references (Wick *et al.* 2025, Autocycler paper) | yes |
| Ion Torrent | `iontorrent_polish_benchmark.sh` | *S. goeteborgense* CCUG 66741ᵀ, *Vibrio* sp. ED-400 | real ONT + Ion PGM (SE) + Illumina, one DNA source each | same-isolate ONT+Illumina assembly | yes |

## Gold standards and their caveats

- **AFsim** — the only track where the truth is mathematically exact; but Medaka on *simulated* reads is
  unrepresentative (badread's error model is not a basecaller's), so judge Medaka on the real tracks.
- **AF2122** — real reads from the strain the reference was finished from: exact truth, rough CLR base
  assembly. A deliberate stress test.
- **MG1655** — reads are from DSM 18039, truth is U00096.3: a handful of real lab-stock differences form an
  error floor no polisher can cross (~2 mismatches + 6 indels). The HiFi arm shows what polishing with
  reads *from a different culture* does to a near-perfect assembly.
- **MBWGS426** — field strain vs the AF2122/97 reference: strain divergence dominates, so only the
  *differences between arms* are meaningful, never the absolute numbers.
- **Non-model** (Autocycler isolates) — same-isolate truth for species far from the usual training
  organisms (*Listeria innocua*: Bacillota; *Providencia rettgeri*: Morganellaceae), chosen to test whether
  Medaka's bacterial-methylation model generalizes beyond model organisms. Caveat: the references were
  produced with Trycycler + Medaka + Polypolish + pyPOLCA (Wick *et al.* 2025), so the short-read polishers
  under test also shaped the truth; the Illumina polishing step scrubs most Medaka-specific errors out of
  the truth, but a shared blind spot cannot be excluded. Accessions (BioProject of the Autocycler paper):
  - *L. innocua* AUSMDU00097349 — BioSample SAMN46906078, ONT SRR32486273, Illumina SRR32484996
  - *P. rettgeri* AUSMDU00026122 — BioSample SAMN46905808, ONT SRR32486258, Illumina SRR32458788
  - References: figshare article 28910099 ("Autocycler: assemblies and reference genomes")

- **Ion Torrent** — does Ion Torrent polishing of a Nanopore assembly fail on exactly the homopolymer
  indels Illumina fixes (shared flow-chemistry error mode)? Three platforms per isolate from one DNA
  submission; the Illumina polish arm is depth-matched to the Ion arm, and every remaining indel is
  classified by homopolymer context. Caveat: the truth is Illumina-polished, so the Ion-vs-none delta
  and the homopolymer classification are the meaningful readouts, not Illumina-vs-truth. Accessions:
  - *S. goeteborgense* CCUG 66741 — SAMN05150228: ONT SRR11858646, Ion SRR8561413, NovaSeq SRR11858647
  - *Vibrio* sp. ED-400 — SAMN20877004: ONT SRR16167661, Ion SRR16167662, MiSeq SRR16167664

## Reproducing

Every script takes `<root_folder> [threads]`, downloads/stages what it can fetch publicly under
`<root_folder>/data_*`, and writes runs + `benchmark_results.tsv` under `<root_folder>/runs_*`. Downloads
go through ENA fastq mirrors (no sra-tools needed) except where a script's header says otherwise.

```bash
conda activate HashBrown
bash validation/benchmarks/autocycler_nonmodel_benchmark.sh /scratch/hashbrown_benchmark 32
```

Results of all tracks, and the defaults that follow from them, are summarized in the wiki's Validation page.
