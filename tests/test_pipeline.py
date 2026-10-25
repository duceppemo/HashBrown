"""End-to-end runs of the three pipelines with stubbed external programs."""

from hashbrown import cli


def run(argv):
    cli.main(argv)


def test_hybrid_pipeline(stub_bin, fastq_gz, tmp_path):
    long_dir, short_dir, out = tmp_path / 'long', tmp_path / 'short', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')
    fastq_gz(long_dir / 'unclassified.fastq.gz')  # must be dropped
    fastq_gz(short_dir / 'sample1_R1.fastq.gz')
    fastq_gz(short_dir / 'sample1_R2.fastq.gz')

    run(['-l', str(long_dir), '-s', str(short_dir), '-o', str(out),
         '--trim-short', '-t', '2', '-p', '1'])

    # Assembly, graphs and stats
    assert (out / '3_assembled_long' / 'sample1.fasta').is_file()
    assert not (out / '3_assembled_long' / 'unclassified.fasta').exists()
    qc = out / '7_assembly_qc'
    assert (qc / 'assembly_graphs' / 'sample1.gfa').is_file()
    assert (qc / 'assembly_graphs' / 'sample1.png').is_file()
    stats = (out / 'flye_stats.tsv').read_text()
    assert 'sample1' in stats and '12000' in stats

    # Polishing: medaka, then polypolish + pypolca with renamed contigs
    assert (out / '4_polished_long' / 'sample1.fasta').is_file()
    final = out / '6_polished_short' / 'sample1.fasta'
    assert final.read_text().startswith('>contig_1\n')
    assert (out / '6_polished_short' / 'sample1.report').is_file()
    assert (out / '5a_trimmed_short' / 'sample1_R1.fastq.gz').is_file()
    assert (out / '5b_trimmed_short_report' / 'sample1.html').is_file()

    # QC
    assert (qc / 'quast' / 'report.tsv').is_file()
    assert (qc / 'mapping_qc_long_reads' / 'sample1.flagstat.txt').is_file()
    assert (qc / 'mapping_qc_long_reads' / 'sample1.mosdepth.summary.txt').is_file()
    assert (qc / 'mapping_qc_short_reads' / 'sample1.flagstat.txt').is_file()
    assert (qc / 'mummer_long' / 'sample1_medaka.ps').is_file()
    assert (qc / 'mummer_short' / 'sample1_polypolish.ps').is_file()
    assert 'sample1\t10' in (qc / 'mapping_qc_long_reads' / 'long_read_cov.tsv').read_text()
    assert 'sample1\t10' in (qc / 'mapping_qc_short_reads' / 'short_read_cov.tsv').read_text()

    # No working leftovers
    assert not (out / '3_assembled_long' / 'sample1').exists()  # flye working folder
    assert not list((out / '6_polished_short').glob('*.sam'))
    assert not list(qc.rglob('*.bam'))
    assert (out / 'log.txt').read_text().startswith('Python')

    # Resuming skips the done steps and does not duplicate anything
    (out / '3_assembled_long' / 'sample1.fasta').write_text('>contig_1\nACGT\n')  # marker survives a resume
    run(['-l', str(long_dir), '-s', str(short_dir), '-o', str(out),
         '--trim-short', '-t', '2', '-p', '1'])
    assert (out / '3_assembled_long' / 'sample1.fasta').read_text() == '>contig_1\nACGT\n'
    assert (qc / 'mapping_qc_long_reads' / 'long_read_cov.tsv').read_text().count('sample1') == 1
    assert 'sample1' in (out / 'flye_stats.tsv').read_text()  # the stats table was not clobbered


def test_long_only_shasta_with_trimming_and_filtering(stub_bin, fastq_gz, tmp_path):
    long_dir, out = tmp_path / 'long', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')

    run(['-l', str(long_dir), '-o', str(out), '-a', 'shasta',
         '--trim-long', '--filter-long', '--min-size', '2000', '-p', '1'])

    assert (out / '1_trimmed_long' / 'sample1.fastq.gz').is_file()
    assert (out / '2_filtered_long' / 'sample1.fastq.gz').is_file()
    assert (out / '3_assembled_long' / 'sample1.fasta').is_file()
    assert (out / '4_polished_long' / 'sample1.fasta').is_file()  # medaka by default
    stats = (out / 'shasta_stats.tsv').read_text()
    assert 'sample1' in stats and '9000000' in stats
    assert (out / '7_assembly_qc' / 'quast' / 'report.tsv').is_file()
    assert not (out / '6_polished_short').exists()  # no short reads


def test_long_only_without_medaka(stub_bin, fastq_gz, tmp_path):
    long_dir, out = tmp_path / 'long', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')

    run(['-l', str(long_dir), '-o', str(out), '--long-polisher', 'none', '-p', '1'])

    assert (out / '3_assembled_long' / 'sample1.fasta').is_file()
    assert not (out / '4_polished_long').exists()
    assert not (out / '7_assembly_qc' / 'mummer_long').exists()  # nothing to compare


def test_hybrid_without_short_polishing(stub_bin, fastq_gz, tmp_path):
    long_dir, short_dir, out = tmp_path / 'long', tmp_path / 'short', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')
    fastq_gz(short_dir / 'sample1_R1.fastq.gz')
    fastq_gz(short_dir / 'sample1_R2.fastq.gz')

    run(['-l', str(long_dir), '-s', str(short_dir), '-o', str(out),
         '--short-polisher', 'none', '-p', '1'])

    assert (out / '4_polished_long' / 'sample1.fasta').is_file()
    assert not (out / '6_polished_short').exists()


def test_short_only_skesa(stub_bin, fastq_gz, tmp_path):
    short_dir, out = tmp_path / 'short', tmp_path / 'out'
    fastq_gz(short_dir / 'pe_R1.fastq.gz')
    fastq_gz(short_dir / 'pe_R2.fastq.gz')
    fastq_gz(short_dir / 'se.fastq.gz')  # single-end sample

    run(['-s', str(short_dir), '-o', str(out), '--trim-short', '-p', '1'])

    assert (out / '1_trimmed_short' / 'pe_R1.fastq.gz').is_file()
    assert (out / '1_trimmed_short' / 'se.fastq.gz').is_file()
    assert (out / '2_assembled_short' / 'pe.fasta').is_file()
    assert (out / '2_assembled_short' / 'se.fasta').is_file()
    qc = out / '3_assembly_qc'
    assert (qc / 'assembly_graphs' / 'pe.gfa').is_file()  # gfa_connector
    assert (qc / 'assembly_graphs' / 'pe.png').is_file()  # Bandage
    assert (qc / 'quast' / 'report.tsv').is_file()
    assert (qc / 'mapping_qc_short_reads' / 'pe.flagstat.txt').is_file()
    assert (qc / 'mapping_qc_short_reads' / 'se.mosdepth.summary.txt').is_file()
    coverage = (qc / 'mapping_qc_short_reads' / 'short_read_cov.tsv').read_text()
    assert 'pe\t10' in coverage and 'se\t10' in coverage


def test_short_only_spades(stub_bin, fastq_gz, tmp_path):
    short_dir, out = tmp_path / 'short', tmp_path / 'out'
    fastq_gz(short_dir / 'pe_R1.fastq.gz')
    fastq_gz(short_dir / 'pe_R2.fastq.gz')

    run(['-s', str(short_dir), '-o', str(out), '-a', 'spades', '--ion-torrent', '-p', '1'])

    assert (out / '2_assembled_short' / 'pe.fasta').is_file()
    assert (out / '3_assembly_qc' / 'assembly_graphs' / 'pe.gfa').is_file()
    assert not (out / '2_assembled_short' / 'pe').exists()  # SPAdes working folder removed
