# Installation

HashBrown is a Python package (≥ 3.10) that drives command-line bioinformatics programs. The package needs
only `psutil`; the programs come from [bioconda](https://bioconda.github.io/).

## With conda (recommended)

```bash
git clone https://github.com/duceppemo/HashBrown
cd HashBrown
conda env create -f environment.yml    # or: mamba env create -f environment.yml
conda activate HashBrown
pip install .

hashbrown --version
```

`environment.yml` installs every program the three pipelines use. If you only ever assemble short reads,
you can trim it down to `skesa`/`spades`, `fastp`, `samtools`, `minimap2`, `quast`, `mosdepth` and
`bandage`.

## Without installing

From a clone with the environment active, `python hashbrown.py …` behaves exactly like the `hashbrown`
command.

## Optional programs

- **Bandage** draws PNG pictures of the assembly graphs.
- **gfa_connector** (shipped with the bioconda SKESA package) rebuilds an assembly graph for SKESA
  assemblies, which have none natively.

If either is missing, HashBrown says so once and skips that output; nothing else changes.

## Troubleshooting

- *The conda solver takes forever or fails.* Use [mamba](https://mamba.readthedocs.io/), and make sure the
  channels are exactly `conda-forge` then `bioconda` (as in `environment.yml`). Historically Medaka's
  TensorFlow builds were the culprit — the reason BGA and BGA_short were separate tools — but Medaka ≥ 2.0
  is PyTorch-based and coexists with SKESA/SPAdes.
- *`Some required programs are not installed or not on the PATH`.* Activate the environment
  (`conda activate HashBrown`). The message lists exactly which programs are missing; HashBrown only
  requires the ones your flags use.
- *Medaka polishing never finishes or crashes.* Check that `bcftools` and `samtools` load
  (`samtools version`): a broken conda install of those is the usual cause. Reinstalling the environment
  fixes it.
