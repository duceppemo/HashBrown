# Contributing

Thanks for helping improve HashBrown. Bug reports, benchmark results and documentation fixes are as welcome
as code.

## Reporting a bug

Open an [issue](https://github.com/duceppemo/HashBrown/issues/new) with the exact command, the `log.txt`
from the output folder (it records the version of every program HashBrown ran) and what the failing step
printed.

## Changing code

```bash
git clone https://github.com/duceppemo/HashBrown && cd HashBrown
conda env create -f environment.yml && conda activate HashBrown
pip install -e ".[test]"
pre-commit install
```

- Keep `ruff check .` and `pytest` green; add a test for every fix or feature. The tests replace the
  external programs with stubs (`tests/conftest.py`), so they run in seconds without conda.
- HashBrown itself needs only `psutil` beyond the Python standard library; external programs are called on
  the command line.
- A new or replacement program (an assembler or a polisher) must be maintained, installable from bioconda
  together with the others, and benchmarked: run `validation/benchmark_polishers.sh` against a gold
  standard and compare the error rates with the current defaults.
- User-visible changes go in `CHANGELOG.md`; options and outputs are documented in `docs/wiki/`.
- Open the pull request against `main`.
