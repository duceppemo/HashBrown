import gzip

import pytest

from hashbrown.methods import Methods, Sample


def test_get_sample_name():
    assert Methods.get_sample_name('sample1_R1.fastq.gz') == 'sample1'
    assert Methods.get_sample_name('/path/to/sample2.fastq.gz') == 'sample2'
    assert Methods.get_sample_name('sample3_barcode01_pass.fq') == 'sample3'
    assert Methods.get_sample_name('plain.fasta') == 'plain'


def test_check_input_folder(tmp_path):
    with pytest.raises(Exception, match='not an existing folder'):
        Methods.check_input_folder(str(tmp_path / 'nope'))

    empty = tmp_path / 'empty'
    empty.mkdir()
    with pytest.raises(Exception, match='No sequence file'):
        Methods.check_input_folder(str(empty))

    (empty / 'reads.fastq.gz').write_bytes(b'')
    Methods.check_input_folder(str(empty))  # no exception

    nested = tmp_path / 'nested'
    (nested / 'sub').mkdir(parents=True)
    (nested / 'sub' / 'reads.fq').write_text('')
    Methods.check_input_folder(str(nested))  # found recursively


def test_check_cpus_clamps_invalid_values():
    from multiprocessing import cpu_count
    total = cpu_count()

    assert Methods.check_cpus(total, 1) == (total, 1)
    assert Methods.check_cpus(0, 2)[0] == total  # 0 threads: clamped to the total
    assert Methods.check_cpus(total * 10, 2)[0] == total  # too many: clamped
    cpu, parallel = Methods.check_cpus(2, 100)  # more parallel samples than threads
    assert parallel <= cpu


def test_check_mem_clamps_invalid_values():
    assert Methods.check_mem(1) == 1
    assert Methods.check_mem(0) >= 1  # 0: set to the default
    assert Methods.check_mem(10 ** 9) < 10 ** 9  # more than the machine has: clamped


def test_get_nanopore_files(tmp_path, capsys):
    folder = tmp_path / 'reads'
    folder.mkdir()
    (folder / 'sample1.fastq.gz').write_bytes(b'')
    (folder / 'sample2_pass.fq').write_text('')
    (folder / 'notes.txt').write_text('ignored')

    sample_dict = Methods.get_nanopore_files(str(folder))
    assert sorted(sample_dict) == ['sample1', 'sample2']
    assert sample_dict['sample1'].nanopore.raw.endswith('sample1.fastq.gz')

    # A second file for the same sample is reported and ignored
    (folder / 'sample1_bis.fastq').write_text('')
    sample_dict = Methods.get_nanopore_files(str(folder))
    assert sample_dict['sample1'].nanopore.raw.endswith('sample1.fastq.gz')
    assert 'More than one long-read file' in capsys.readouterr().err

    empty = tmp_path / 'no_reads'
    empty.mkdir()
    with pytest.raises(Exception, match='No sequence file'):
        Methods.get_nanopore_files(str(empty))


def test_get_illumina_files_hybrid(tmp_path, capsys):
    long_folder = tmp_path / 'long'
    long_folder.mkdir()
    (long_folder / 'sample1.fastq.gz').write_bytes(b'')
    short_folder = tmp_path / 'short'
    short_folder.mkdir()
    (short_folder / 'sample1_R1.fastq.gz').write_bytes(b'')
    (short_folder / 'sample1_R2.fastq.gz').write_bytes(b'')
    (short_folder / 'orphan_R1.fastq.gz').write_bytes(b'')  # no matching long reads
    (short_folder / 'sample1.fastq.gz').write_bytes(b'')  # a stray third file: the pair wins

    sample_dict = Methods.get_nanopore_files(str(long_folder))
    sample_dict = Methods.get_illumina_files(str(short_folder), sample_dict)

    assert sorted(sample_dict) == ['sample1']  # the orphan did not create a sample
    r1, r2 = sample_dict['sample1'].illumina.raw
    assert '_R1' in r1 and '_R2' in r2
    err = capsys.readouterr().err
    assert 'No matching long reads' in err
    assert 'Too many short-read files' in err


def test_get_illumina_files_hybrid_single_end(tmp_path):
    """A single-end short-read file (e.g. Ion Torrent) is attached to its long-read sample."""
    long_folder = tmp_path / 'long'
    long_folder.mkdir()
    (long_folder / 'sample1.fastq.gz').write_bytes(b'')
    short_folder = tmp_path / 'short'
    short_folder.mkdir()
    (short_folder / 'sample1.fastq.gz').write_bytes(b'')

    sample_dict = Methods.get_nanopore_files(str(long_folder))
    sample_dict = Methods.get_illumina_files(str(short_folder), sample_dict)

    assert len(sample_dict['sample1'].illumina.raw) == 1


def test_get_illumina_files_standalone(tmp_path):
    short_folder = tmp_path / 'short'
    short_folder.mkdir()
    (short_folder / 'pe_R1.fastq.gz').write_bytes(b'')
    (short_folder / 'pe_R2.fastq.gz').write_bytes(b'')
    (short_folder / 'se.fastq.gz').write_bytes(b'')

    sample_dict = Methods.get_illumina_files(str(short_folder), dict(), require_long_reads=False)
    assert sorted(sample_dict) == ['pe', 'se']
    assert len(sample_dict['pe'].illumina.raw) == 2
    assert len(sample_dict['se'].illumina.raw) == 1


def test_sample_best_helpers():
    sample = Sample()
    sample.nanopore.raw = 'raw.fq'
    assert sample.best_nanopore() == 'raw.fq'
    sample.nanopore.trimmed = 'trimmed.fq'
    assert sample.best_nanopore() == 'trimmed.fq'
    sample.nanopore.filtered = 'filtered.fq'
    assert sample.best_nanopore() == 'filtered.fq'

    assert sample.best_illumina() == []
    sample.illumina.raw = ['r1', 'r2']
    assert sample.best_illumina() == ['r1', 'r2']
    sample.illumina.trimmed = ['t1', 't2']
    assert sample.best_illumina() == ['t1', 't2']

    assert sample.best_assembly() == ''
    sample.assembly.raw = 'raw.fasta'
    sample.assembly.medaka = 'medaka.fasta'
    sample.assembly.polypolish = 'final.fasta'
    assert sample.best_assembly() == 'final.fasta'


def test_format_fasta(tmp_path):
    input_fasta = tmp_path / 'in.fasta'
    output_fasta = tmp_path / 'out.fasta'
    input_fasta.write_text('>contig_1 info\n' + 'A' * 100 + '\n' + 'C' * 10 + '\n\n>contig_2\nACGT\n')

    Methods.format_fasta(str(input_fasta), str(output_fasta))

    lines = output_fasta.read_text().splitlines()
    assert lines[0] == '>contig_1 info'
    assert lines[1] == 'A' * 80
    assert lines[2] == 'A' * 20 + 'C' * 10
    assert lines[3] == '>contig_2'
    assert lines[4] == 'ACGT'


def test_list_to_file_and_flag_done(tmp_path):
    out = tmp_path / 'list.txt'
    Methods.list_to_file(['a', 'b'], str(out))
    assert out.read_text() == 'a\nb\n'

    flag = tmp_path / 'done'
    Methods.flag_done(str(flag))
    assert flag.exists()


def test_fasta_length(tmp_path):
    from hashbrown.assembly_qc import AssemblyQcMethods

    assert AssemblyQcMethods.fasta_length('') == 0
    assert AssemblyQcMethods.fasta_length(str(tmp_path / 'nope.fasta')) == 0

    plain = tmp_path / 'a.fasta'
    plain.write_text('>c1\nACGT\nACGT\n>c2\nAC\n')
    assert AssemblyQcMethods.fasta_length(str(plain)) == 10

    gz = tmp_path / 'a.fasta.gz'
    with gzip.open(gz, 'wt') as f:
        f.write('>c1\nACGTACGT\n')
    assert AssemblyQcMethods.fasta_length(str(gz)) == 8
