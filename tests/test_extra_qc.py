import os

from hashbrown.extra_qc import ExtraQcMethods
from hashbrown.methods import Sample


def test_modules_note_missing_tools(tmp_path, monkeypatch):
    """Without the optional tools installed, every module returns empty results and a hint."""
    monkeypatch.setenv('PATH', str(tmp_path))  # nothing on the PATH

    results, notes = ExtraQcMethods.species_ani('s', 'asm.fasta', str(tmp_path), 1)
    assert results == {} and 'mash' in notes[0]

    results, notes = ExtraQcMethods.checkm2('s', 'asm.fasta', str(tmp_path), 1)
    assert results == {} and 'CheckM2' in notes[0]

    results, notes = ExtraQcMethods.busco('s', 'asm.fasta', str(tmp_path), 1)
    assert results == {} and 'BUSCO' in notes[0]

    results, notes = ExtraQcMethods.rrna_trna('s', 'asm.fasta', str(tmp_path), 1)
    assert results == {} and len(notes) == 2


def test_merqury_skipped_without_short_reads(tmp_path):
    qv, notes = ExtraQcMethods.merqury_qv('s', Sample(), [], str(tmp_path), 1)
    assert qv == {} and notes == []


def test_species_ani_without_sketch(stub_bin, tmp_path, monkeypatch):
    monkeypatch.setenv('HASHBROWN_DB', str(tmp_path / 'empty_db'))
    results, notes = ExtraQcMethods.species_ani('s', 'asm.fasta', str(tmp_path), 1)
    assert results == {} and 'hashbrown db mash' in notes[0]


def test_phase3_modules_note_missing_tools(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path))

    results, notes = ExtraQcMethods.mlst('s', 'asm.fasta')
    assert results == {} and 'mlst' in notes[0]

    rows, notes = ExtraQcMethods.mob_suite('s', 'asm.fasta', str(tmp_path), 1)
    assert rows == [] and 'MOB-suite' in notes[0]

    rows, notes = ExtraQcMethods.amr_virulence('s', 'asm.fasta', str(tmp_path), 1)
    assert rows == [] and len(notes) == 2  # AMRFinderPlus and abricate

    results, notes = ExtraQcMethods.annotate('s', 'asm.fasta', str(tmp_path), 1, mode='auto')
    assert results == {} and 'Bakta' in notes[0]

    results, notes = ExtraQcMethods.annotate('s', 'asm.fasta', str(tmp_path), 1, mode='none')
    assert results == {} and notes == []

    results, notes = ExtraQcMethods.annotate('s', 'asm.fasta', str(tmp_path), 1, mode='pgap',
                                             organism='Escherichia coli')
    assert results == {} and 'pgap.py' in notes[0]


def test_bakta_needs_database(stub_bin, tmp_path, monkeypatch):
    monkeypatch.setenv('HASHBROWN_DB', str(tmp_path / 'empty_db'))
    monkeypatch.delenv('BAKTA_DB', raising=False)
    results, notes = ExtraQcMethods.annotate('s', 'asm.fasta', str(tmp_path), 1, mode='bakta')
    assert results == {} and 'hashbrown db bakta' in notes[0]


def test_reorient_skipped_without_dnaapler(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path))
    done_flag = tmp_path / 'done_reorienting'
    ExtraQcMethods.reorient_parallel({'s': Sample()}, str(done_flag), 1)
    assert not done_flag.exists()  # so it runs once dnaapler is installed


def test_reorient_in_place_keeps_sidecar(stub_bin, tmp_path):
    fasta = tmp_path / 'sample1.fasta'
    fasta.write_text('>contig_1\nACGTACGT\n')
    sample_obj = Sample()
    sample_obj.assembly.raw = str(fasta)

    assert ExtraQcMethods.reorient_dnaapler('sample1', sample_obj, 1)
    assert (tmp_path / 'sample1_unoriented.fasta').read_text() == '>contig_1\nACGTACGT\n'
    assert fasta.read_text().startswith('>contig_1')
    assert sample_obj.assembly.unoriented == str(tmp_path / 'sample1_unoriented.fasta')
    assert not os.path.isdir(tmp_path / 'sample1_dnaapler')  # working folder cleaned up
