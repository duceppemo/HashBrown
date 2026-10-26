import re
from pathlib import Path

import pytest

import hashbrown
from hashbrown import cli

REPO = Path(__file__).resolve().parents[1]


def test_version_in_sync_with_pyproject():
    # tomllib only exists from Python 3.11; the project supports 3.10
    match = re.search(r'^version = "([^"]+)"$', (REPO / 'pyproject.toml').read_text(), re.MULTILINE)
    assert match is not None
    assert match.group(1) == hashbrown.__version__


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.parse_arguments(['--version'])
    assert excinfo.value.code == 0
    assert hashbrown.__version__ in capsys.readouterr().out


@pytest.mark.parametrize('argv, message', [
    (['-l', 'long'], 'output folder'),
    (['-o', 'out'], 'provide long reads'),
    (['-l', 'long', '-o', 'out', '-a', 'skesa'], 'assembles short reads only'),
    (['-l', 'long', '-o', 'out', '-a', 'spades'], 'assembles short reads only'),
    (['-s', 'short', '-o', 'out', '-a', 'flye'], 'needs long reads'),
    (['-s', 'short', '-o', 'out', '-a', 'shasta'], 'needs long reads'),
    (['-s', 'short', '-o', 'out', '--size', '5000000'], 'need long reads'),
    (['-s', 'short', '-o', 'out', '--ion-torrent'], 'need SPAdes'),
])
def test_incoherent_flags_are_refused(argv, message, capsys):
    with pytest.raises(SystemExit):
        cli.parse_arguments(argv)
    assert message in capsys.readouterr().err


def test_default_assembler_follows_the_inputs():
    assert cli.parse_arguments(['-l', 'long', '-o', 'out']).assembler == 'flye'
    assert cli.parse_arguments(['-s', 'short', '-o', 'out']).assembler == 'skesa'
    assert cli.parse_arguments(['-l', 'long', '-o', 'out', '-a', 'shasta']).assembler == 'shasta'
    assert cli.parse_arguments(['-s', 'short', '-o', 'out', '-a', 'spades', '--ion-torrent']).ion_torrent
    # Ion Torrent is also valid in hybrid mode (the reads are preprocessed differently and a
    # polishing warning is printed at run time)
    assert cli.parse_arguments(['-l', 'long', '-s', 'short', '-o', 'out', '--ion-torrent']).ion_torrent


def test_read_prep_defaults():
    args = cli.parse_arguments(['-l', 'long', '-s', 'short', '-o', 'out'])
    assert not args.trim_long  # Porechop is opt-in: unmaintained, slow, redundant on recent basecalls
    assert args.filter_long and args.trim_short
    args = cli.parse_arguments(['-l', 'long', '-o', 'out', '--trim-long', '--no-filter-long'])
    assert args.trim_long and not args.filter_long
    # short-only runs ignore the long-read prep steps instead of erroring
    args = cli.parse_arguments(['-s', 'short', '-o', 'out'])
    assert not args.trim_long and not args.filter_long and args.trim_short


def test_short_polisher_choices():
    args = cli.parse_arguments(['-l', 'long', '-o', 'out'])
    assert cli.SHORT_POLISHERS[args.short_polisher] == ('polypolish', 'pypolca')
    args = cli.parse_arguments(['-l', 'long', '-o', 'out', '--short-polisher', 'none'])
    assert cli.SHORT_POLISHERS[args.short_polisher] == ()


def test_list_models(stub_bin, capsys):
    cli.main(['--list-models'])
    out = capsys.readouterr().out
    assert 'model_a' in out and 'r1041_e82_400bps_sup_g615' in out


def test_list_models_without_medaka(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('PATH', str(tmp_path))  # nothing installed
    cli.main(['--list-models'])
    assert 'Could not list' in capsys.readouterr().out


def test_missing_input_folder_is_a_clean_error(tmp_path, stub_bin):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['-l', str(tmp_path / 'nope'), '-o', str(tmp_path / 'out')])
    assert 'not an existing folder' in str(excinfo.value)


def test_missing_program_is_a_clean_error(tmp_path, monkeypatch, fastq_gz):
    fastq_gz(tmp_path / 'long' / 'sample1.fastq.gz')
    monkeypatch.setenv('PATH', str(tmp_path))  # nothing installed
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['-l', str(tmp_path / 'long'), '-o', str(tmp_path / 'out')])
    assert 'not installed or not on the PATH' in str(excinfo.value)
