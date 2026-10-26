from hashbrown.report import ReportMethods


def test_combine_stats_weights_rates_by_bases():
    rows = [{'num_seqs': '10', 'sum_len': '1000', 'N50': '100', 'AvgQual': '30.0', 'GC(%)': '40.0'},
            {'num_seqs': '10', 'sum_len': '3000', 'N50': '300', 'AvgQual': '34.0', 'GC(%)': '60.0'}]
    combined = ReportMethods._combine_stats(rows)
    assert combined['reads'] == 20
    assert combined['bases'] == 4000
    assert combined['n50'] == 300
    assert combined['avg_qual'] == 33.0  # (30*1000 + 34*3000) / 4000
    assert combined['gc'] == 55.0


def test_combine_stats_empty():
    assert ReportMethods._combine_stats([]) is None


def test_parse_mosdepth_summary(tmp_path):
    summary = tmp_path / 's.mosdepth.summary.txt'
    summary.write_text('chrom\tlength\tbases\tmean\tmin\tmax\n'
                       'chr1\t1000\t50000\t50.00\t1\t80\n'
                       'chr1_region\t1000\t50000\t50.00\t1\t80\n'
                       'plasmid\t100\t2000\t20.00\t1\t40\n'
                       'total\t1100\t52000\t47.27\t1\t80\n')
    contigs, total = ReportMethods.parse_mosdepth_summary(str(summary))
    assert [c['contig'] for c in contigs] == ['chr1', 'plasmid']
    assert contigs[0]['mean'] == 50.0
    assert total['length'] == 1100


def test_parse_mosdepth_summary_missing():
    contigs, total = ReportMethods.parse_mosdepth_summary('/nonexistent')
    assert contigs == [] and total is None


def test_sample_reads(fastq_gz, tmp_path):
    fastq = fastq_gz(tmp_path / 'reads.fastq.gz', n_reads=5)
    lengths, quals = ReportMethods.sample_reads(str(fastq))
    assert lengths == [12] * 5
    assert quals == [40.0] * 5  # 'I' is Phred 40


def test_fig_read_distributions_empty():
    assert ReportMethods.fig_read_distributions([('raw', [], [])]) == ''
