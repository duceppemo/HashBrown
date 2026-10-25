# Development

## Layout

```
hashbrown/            the package
├── cli.py            argument parsing + the three pipelines (HashBrown.run_long / run_short)
├── methods.py        input discovery, Sample, resource checks, small file utilities
├── nanopore.py       Porechop, Filtlong, Flye, Shasta, Medaka
├── illumina.py       fastp, SKESA, SPAdes, Polypolish, pyPOLCA
└── assembly_qc.py    minimap2 mapping, Qualimap, QUAST, nucmer dot plots, coverage
tests/                pytest suite; conftest.py stubs every external program
validation/           polisher benchmark against a gold standard
docs/wiki/            these pages (published by .github/workflows/wiki.yml)
docs/images/          logo and pipeline diagram (regenerate with make_svgs.py)
recipe/               bioconda recipe
```

## Tests

```bash
pip install -e ".[test]"
pytest
ruff check .
```

The tests run in seconds without conda: `tests/conftest.py` builds a folder of fake `flye`, `medaka_consensus`,
`samtools`, … executables that mimic just enough of each program's command line, and the pipeline tests run
the real orchestration end to end against them (including resuming). A change to how a program is *called*
must update its stub to match the real program's behaviour.

## CI and releases

- `ci.yml`: ruff + pytest on Linux (Python 3.10–3.13) and macOS, coverage to Codecov (needs the
  `CODECOV_TOKEN` secret).
- `wiki.yml`: pushes `docs/wiki/` to the GitHub wiki on every change (create one page in the web UI once,
  so the wiki repository exists).
- `release.yml`: on a `v*` tag — checks the tag matches `hashbrown.__version__` and `pyproject.toml`,
  builds the wheel/sdist, and creates a GitHub release with the matching `CHANGELOG.md` section.

Releasing: bump the version in `hashbrown/__init__.py` **and** `pyproject.toml`, date the section in
`CHANGELOG.md`, update `CITATION.cff`, then `git tag v<version> && git push --tags`. For bioconda, update
`recipe/meta.yaml` with the new tarball sha256 and submit to bioconda-recipes.

## Adding an assembler or polisher

1. Add the runner in `illumina.py`/`nanopore.py` and wire it in `cli.py` (choices + dependency map in
   `methods.py`).
2. Add its stub in `tests/conftest.py` and a pipeline test.
3. Benchmark it (`validation/benchmark_polishers.sh`) — a polisher that adds errors in low-coverage
   regions does not go in, no matter its reputation.
4. Document it in `docs/wiki/Methods.md` and `CHANGELOG.md`, and add it to `environment.yml` and
   `recipe/meta.yaml`.
