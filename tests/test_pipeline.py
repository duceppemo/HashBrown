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
    assert 'sample1\traw\t10' in (qc / 'mapping_qc_long_reads' / 'long_read_cov.tsv').read_text()
    short_cov = (qc / 'mapping_qc_short_reads' / 'short_read_cov.tsv').read_text()
    assert 'sample1\traw\t10' in short_cov and 'sample1\ttrimmed\t10' in short_cov  # one row per stage
    assert (qc / 'mapping_qc_short_reads' / 'sample1_trimmed.flagstat.txt').is_file()
    assert (qc / 'mapping_qc_short_reads' / 'sample1_raw.mosdepth.summary.txt').is_file()

    # Report: per-sample page, run index and the machine-readable tables
    report = out / '8_report'
    page = (report / 'sample1.html').read_text()
    assert 'Read QC - Long reads' in page and 'Read QC - Short reads' in page
    assert 'Assembly - stage by stage' in page and 'Polished - short reads' in page
    assert 'Coverage per contig' in page
    assert 'sample1.html' in (report / 'index.html').read_text()
    funnel = (report / 'tables' / 'sample1_read_funnel_long_reads.tsv').read_text()
    assert funnel.startswith('stage\treads\tbases')
    assert (report / 'tables' / 'sample1_assembly_funnel.tsv').is_file()
    assert (report / 'tables' / 'sample1_contig_coverage.tsv').is_file()

    # Optional QC modules all ran through their stubs
    assert (out / '6_polished_short' / 'sample1_unoriented.fasta').is_file()  # dnaapler, in place
    assert 'Reoriented (dnaapler)' in page
    assert '54.50' in page  # Merqury QV column
    assert 'Species confirmation' in page and 'Escherichia coli' in page and '99.99%' in page
    assert 'CheckM2 completeness' in page and '99.98' in page
    assert 'bacteria_odb10' in page  # BUSCO
    assert '1x 16S' in page and '2 tRNAs' in page
    summary = (report / 'tables' / 'sample1_identity_completeness.tsv').read_text()
    assert 'completeness' in summary and 'ani' in summary

    # Phase 3: typing, genes, annotation
    assert 'ecoli' in page and 'ST' in page  # MLST, also in the index
    assert 'ecoli ST10' in (report / 'index.html').read_text()
    assert 'plasmid_AA474' in page and 'IncFIB' in page  # MOB-suite
    assert 'blaEC' in page and 'gene_vfdb' in page and 'gene_plasmidfinder' in page
    assert 'Bakta: 4300 CDSs, 86 tRNAs, 22 rRNAs' in page
    genes = (report / 'tables' / 'sample1_amr_virulence.tsv').read_text()
    assert 'AMRFinderPlus' in genes and 'VFDB' in genes and 'PlasmidFinder' in genes
    assert (report / 'tables' / 'sample1_plasmids.tsv').is_file()

    # Phase 4: GenomeScope2 profile, ConFindr verdict and the blob plot
    assert 'GenomeScope2' in page and '4,500,000 bp - 4,600,000 bp' in page
    assert 'ConFindr verdict' in page and 'Clean' in page
    assert 'Blob plot' in page
    blob = (report / 'tables' / 'sample1_blob.tsv').read_text()
    assert 'contig_1' in blob and 'Escherichia coli' in blob
    assert 'confindr_status' in summary or 'confindr_status' in \
        (report / 'tables' / 'sample1_identity_completeness.tsv').read_text()

    # QA additions: criteria table with status chips, MIMAG verdict, traceability appendix
    assert 'QC criteria' in page and 'status-pass' in page
    assert 'Medium-quality draft (MIMAG)' in page  # the stub's 2 tRNA isotypes keep it below high
    assert 'Traceability' in page and 'Software versions' in page and 'Command line' in page
    assert 'hashbrown ' in page  # the echoed command line
    assert 'Python' in page  # from log.txt

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
    # one row per stage (raw + the default filtering), not duplicated by the resume
    assert (qc / 'mapping_qc_long_reads' / 'long_read_cov.tsv').read_text().count('sample1') == 2
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
    assert 'pe\traw\t10' in coverage and 'se\traw\t10' in coverage

    # Report: one page per sample, both in the index, no long-read section
    report = out / '4_report'
    page = (report / 'se.html').read_text()
    assert 'Read QC - Short reads' in page and 'Read QC - Long reads' not in page
    index = (report / 'index.html').read_text()
    assert 'pe.html' in index and 'se.html' in index


def test_short_only_spades(stub_bin, fastq_gz, tmp_path):
    short_dir, out = tmp_path / 'short', tmp_path / 'out'
    fastq_gz(short_dir / 'pe_R1.fastq.gz')
    fastq_gz(short_dir / 'pe_R2.fastq.gz')

    run(['-s', str(short_dir), '-o', str(out), '-a', 'spades', '--ion-torrent', '-p', '1'])

    assert (out / '2_assembled_short' / 'pe.fasta').is_file()
    assert (out / '3_assembly_qc' / 'assembly_graphs' / 'pe.gfa').is_file()
    assert not (out / '2_assembled_short' / 'pe').exists()  # SPAdes working folder removed


def test_hybrid_ion_torrent_warns_about_polishing(stub_bin, fastq_gz, tmp_path, capsys):
    long_dir, short_dir, out = tmp_path / 'long', tmp_path / 'short', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')
    fastq_gz(short_dir / 'sample1_R1.fastq.gz')
    fastq_gz(short_dir / 'sample1_R2.fastq.gz')

    run(['-l', str(long_dir), '-s', str(short_dir), '-o', str(out), '--ion-torrent', '-p', '1'])

    output = capsys.readouterr().out
    assert 'polishing Nanopore assemblies with Ion Torrent reads is not recommended' in output
    assert (out / '6_polished_short' / 'sample1.fasta').is_file()  # warned, not blocked


def test_hybrid_single_end_short_reads_polish(stub_bin, fastq_gz, tmp_path):
    """Single-end short reads (e.g. Ion Torrent) can polish a long-read assembly."""
    long_dir, short_dir, out = tmp_path / 'long', tmp_path / 'short', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')
    fastq_gz(short_dir / 'sample1.fastq.gz')  # no _R1/_R2: single-end

    run(['-l', str(long_dir), '-s', str(short_dir), '-o', str(out), '-p', '1'])

    final = out / '6_polished_short' / 'sample1.fasta'
    assert final.is_file() and final.read_text().startswith('>contig_1')
    assert (out / '6_polished_short' / 'sample1.report').is_file()  # pypolca ran single-end


def test_skip_report(stub_bin, fastq_gz, tmp_path):
    long_dir, out = tmp_path / 'long', tmp_path / 'out'
    fastq_gz(long_dir / 'sample1.fastq.gz')
    run(['-l', str(long_dir), '-o', str(out), '--skip-report', '-p', '1'])
    assert not (out / '8_report').exists()
