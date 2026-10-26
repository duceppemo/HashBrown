"""Per-sample HTML/PDF report: the effect of every pipeline stage on the reads and the assembly.

Core tools (seqkit, mosdepth, QUAST) are hard requirements of the pipeline; everything else degrades
gracefully: a missing optional tool (pycoQC-ng, WeasyPrint) is reported as skipped, never an error.
"""

import csv
import gzip
import io
import os
import shutil
import subprocess
from datetime import date

from jinja2 import Template

from hashbrown import __version__
from hashbrown.methods import Methods

SAMPLED_READS = 10000  # reads per stage used for the distribution figures


class ReportMethods:
    @staticmethod
    def seqkit_stats(files, cpu):
        """Run `seqkit stats -aT` on the files; returns one dict per file keyed on the header names."""
        files = [f for f in files if f and os.path.exists(f)]
        if not files:
            return []
        cmd = ['seqkit', 'stats', '-a', '-T', '-j', str(cpu)] + files
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return []
        rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter='\t'))
        return rows

    @staticmethod
    def _combine_stats(rows):
        """Combine per-file seqkit rows (e.g. R1+R2) into one: sums for counts, weighted means for rates."""
        if not rows:
            return None
        reads = sum(int(r['num_seqs']) for r in rows)
        bases = sum(int(r['sum_len']) for r in rows)

        def weighted(key):
            try:
                return sum(float(r[key]) * int(r['sum_len']) for r in rows) / bases if bases else 0.0
            except (KeyError, ValueError):
                return None
        return {'reads': reads, 'bases': bases,
                'n50': max(int(r.get('N50', 0)) for r in rows),
                'avg_qual': weighted('AvgQual'), 'gc': weighted('GC(%)')}

    @staticmethod
    def read_funnel(sample, sample_obj, genome_size, cpu):
        """One table per read type: raw -> trimmed -> filtered, with the effect of each stage."""
        funnels = []
        read_sets = [('Long reads', [('raw', [sample_obj.nanopore.raw]),
                                     ('trimmed (Porechop)', [sample_obj.nanopore.trimmed]),
                                     ('filtered (Filtlong)', [sample_obj.nanopore.filtered])]),
                     ('Short reads', [('raw', list(sample_obj.illumina.raw)),
                                      ('trimmed (fastp)', list(sample_obj.illumina.trimmed))])]
        for title, stages in read_sets:
            rows = []
            raw_bases = None
            for stage_name, files in stages:
                stats = ReportMethods._combine_stats(ReportMethods.seqkit_stats(files, cpu))
                if not stats:
                    continue
                if raw_bases is None:
                    raw_bases = stats['bases']
                rows.append({
                    'stage': stage_name,
                    'reads': stats['reads'],
                    'bases': stats['bases'],
                    'n50': stats['n50'],
                    'avg_qual': f"{stats['avg_qual']:.1f}" if stats['avg_qual'] else '-',
                    'gc': f"{stats['gc']:.1f}" if stats['gc'] else '-',
                    'coverage': f"{stats['bases'] / genome_size:.0f}x" if genome_size else '-',
                    'kept': f"{100 * stats['bases'] / raw_bases:.1f}%" if raw_bases else '-',
                    'files': files,
                })
            if rows:
                funnels.append({'title': title, 'rows': rows})
        return funnels

    @staticmethod
    def assembly_funnel(sample, stage_fastas, cpu):
        """One row per assembly stage that exists on disk: contigs, size, N50, GC."""
        rows = []
        for stage_name, fasta in stage_fastas:
            if not fasta or not os.path.exists(fasta):
                continue
            stats = ReportMethods.seqkit_stats([fasta], cpu)
            if not stats:
                continue
            r = stats[0]
            rows.append({'stage': stage_name,
                         'contigs': int(r['num_seqs']),
                         'size': int(r['sum_len']),
                         'n50': int(r.get('N50', 0)),
                         'gc': r.get('GC(%)', '-'),
                         'qv': '-'})
        return rows

    @staticmethod
    def sample_reads(fastq, max_reads=SAMPLED_READS):
        """Lengths and mean Phred scores of the first reads of a fastq(.gz): enough for a distribution."""
        lengths, quals = [], []
        opener = gzip.open if fastq.endswith('.gz') else open
        try:
            with opener(fastq, 'rt') as f:
                for i, line in enumerate(f):
                    if i // 4 >= max_reads:
                        break
                    if i % 4 == 1:
                        lengths.append(len(line.strip()))
                    elif i % 4 == 3:
                        q = line.strip()
                        if q:
                            quals.append(sum(ord(c) - 33 for c in q) / len(q))
        except OSError:
            return [], []
        return lengths, quals

    @staticmethod
    def _inline_svg(fig):
        """The figure as an SVG string ready to embed in HTML (XML prolog and doctype stripped)."""
        from matplotlib import pyplot as plt
        svg = io.StringIO()
        fig.savefig(svg, format='svg')
        plt.close(fig)
        text = svg.getvalue()
        return text[text.find('<svg'):]

    @staticmethod
    def fig_read_distributions(stage_samples):
        """Overlaid per-stage read length and quality distributions; returns an inline SVG string."""
        stage_samples = [(name, lengths, quals) for name, lengths, quals in stage_samples if lengths]
        if not stage_samples:
            return ''
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt

        fig, (ax_len, ax_q) = plt.subplots(1, 2, figsize=(9, 3.2))
        for name, lengths, quals in stage_samples:
            ax_len.hist(lengths, bins=60, histtype='step', linewidth=1.4, label=name)
            if quals:
                ax_q.hist(quals, bins=40, histtype='step', linewidth=1.4, label=name)
        ax_len.set_xlabel('Read length (bp)')
        ax_len.set_ylabel('Reads')
        longest = max(max(lengths) for _, lengths, _ in stage_samples)
        shortest = min(min(lengths) for _, lengths, _ in stage_samples)
        if shortest > 0 and longest > 20 * shortest:
            ax_len.set_xscale('log')
        ax_q.set_xlabel('Mean read quality (Phred)')
        ax_len.legend(frameon=False, fontsize=8)
        for ax in (ax_len, ax_q):
            ax.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        return ReportMethods._inline_svg(fig)

    @staticmethod
    def parse_mosdepth_summary(summary_file):
        """Per-contig mean coverage rows plus the total, from a <sample>.mosdepth.summary.txt."""
        contigs, total = [], None
        if not summary_file or not os.path.exists(summary_file):
            return contigs, total
        with open(summary_file) as f:
            for fields in csv.reader(f, delimiter='\t'):
                if fields[0] == 'chrom' or fields[0].endswith('_region'):
                    continue
                row = {'contig': fields[0], 'length': int(fields[1]), 'mean': float(fields[3])}
                if fields[0] == 'total':
                    total = row
                else:
                    contigs.append(row)
        return contigs, total

    @staticmethod
    def fig_contig_coverage(contigs_long, contigs_short):
        """Grouped bars of per-contig mean coverage, long vs short reads; inline SVG string."""
        names = [c['contig'] for c in (contigs_long or contigs_short)]
        if not names:
            return ''
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt

        long_map = {c['contig']: c['mean'] for c in contigs_long}
        short_map = {c['contig']: c['mean'] for c in contigs_short}
        x = range(len(names))
        width = 0.38 if (long_map and short_map) else 0.76
        fig, ax = plt.subplots(figsize=(max(4.5, 0.9 * len(names) + 2), 3.0))
        if long_map:
            ax.bar([i - width / 2 for i in x] if short_map else list(x),
                   [long_map.get(n, 0) for n in names], width, label='Long reads')
        if short_map:
            ax.bar([i + width / 2 for i in x] if long_map else list(x),
                   [short_map.get(n, 0) for n in names], width, label='Short reads')
        ax.set_xticks(list(x))
        ax.set_xticklabels([n[:20] for n in names], rotation=30, ha='right', fontsize=8)
        ax.set_ylabel('Mean coverage (x)')
        ax.legend(frameon=False, fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        return ReportMethods._inline_svg(fig)

    @staticmethod
    def pycoqc_available():
        return bool(shutil.which('pycoQC') and shutil.which('Fastq_to_seq_summary'))

    @staticmethod
    def run_pycoqc(sample, tag, fastq, pycoqc_folder, cpu):
        """Per-stage pycoQC-ng report via its fastq converter; returns the html path or ''."""
        if not ReportMethods.pycoqc_available():
            return ''
        Methods.make_folder(pycoqc_folder)
        summary = os.path.join(pycoqc_folder, f'{sample}_{tag}_summary.txt.tmp')
        html = os.path.join(pycoqc_folder, f'{sample}_{tag}.html')
        if os.path.exists(html):
            return html
        # The converter walks a directory, so give it one holding just this file
        fastq_dir = os.path.join(pycoqc_folder, f'.{sample}_{tag}_fastq')
        Methods.make_folder(fastq_dir)
        link = os.path.join(fastq_dir, os.path.basename(fastq))
        if not os.path.exists(link):
            os.symlink(os.path.realpath(fastq), link)
        done = subprocess.run(['Fastq_to_seq_summary', '-f', fastq_dir, '-s', summary, '-t', str(cpu)],
                              stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if done.returncode == 0 and os.path.exists(summary) and os.path.getsize(summary):
            subprocess.run(['pycoQC', '-f', summary, '-o', html],
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        for leftover in (summary, link):
            if os.path.exists(leftover):
                os.remove(leftover)
        if os.path.isdir(fastq_dir):
            os.rmdir(fastq_dir)
        return html if os.path.exists(html) else ''

    @staticmethod
    def write_tables(tables_folder, sample, funnels, assembly_rows, contigs_long, contigs_short,
                     identity=None, completeness=None, sanity=None,
                     typing=None, plasmids=None, genes=None):
        """Every report table as machine-readable TSV."""
        Methods.make_folder(tables_folder)

        def write(name, fieldnames, rows):
            with open(os.path.join(tables_folder, f'{sample}_{name}.tsv'), 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                        extrasaction='ignore', restval='-')
                writer.writeheader()
                writer.writerows(rows)
        for funnel in funnels:
            tag = 'long_reads' if funnel['title'].startswith('Long') else 'short_reads'
            write(f'read_funnel_{tag}', ['stage', 'reads', 'bases', 'n50', 'avg_qual', 'gc',
                                         'coverage', 'kept'], funnel['rows'])
        if assembly_rows:
            write('assembly_funnel', ['stage', 'contigs', 'size', 'n50', 'gc', 'qv'], assembly_rows)
        summary = {**(identity or {}), **(completeness or {}), **(sanity or {}), **(typing or {})}
        if summary:
            write('identity_completeness', list(summary), [summary])
        if plasmids:
            write('plasmids', ['plasmid', 'size', 'rep', 'relaxase', 'mobility'], plasmids)
        if genes:
            write('amr_virulence', ['gene', 'source', 'identity', 'coverage', 'contig', 'info'], genes)
        cov_rows = [{'contig': c['contig'], 'length': c['length'], 'read_type': t, 'mean_coverage': c['mean']}
                    for t, contigs in (('long', contigs_long), ('short', contigs_short)) for c in contigs]
        if cov_rows:
            write('contig_coverage', ['contig', 'length', 'read_type', 'mean_coverage'], cov_rows)

    @staticmethod
    def render_pdf(html_file, notes):
        """WeasyPrint render of the html twin; optional."""
        try:
            from weasyprint import HTML
        except Exception:
            notes.append('PDF skipped: WeasyPrint is not installed (conda install weasyprint).')
            return ''
        pdf_file = html_file[:-5] + '.pdf'
        try:
            HTML(filename=html_file).write_pdf(pdf_file)
        except Exception as error:
            notes.append(f'PDF skipped: WeasyPrint failed ({error}).')
            return ''
        return pdf_file

    @staticmethod
    def build_report(sample, sample_obj, report_folder, stage_fastas, mosdepth_long, mosdepth_short,
                     links, cpu, annotate='auto'):
        """Collect everything for one sample and render <sample>.html (+ .pdf); returns a summary dict."""
        notes = []
        genome_size = 0
        final_assembly = sample_obj.best_assembly()
        if final_assembly and os.path.exists(final_assembly):
            from hashbrown.assembly_qc import AssemblyQcMethods
            genome_size = AssemblyQcMethods.fasta_length(final_assembly)

        # A dnaapler sidecar means the last stage file was reoriented in place: show both
        if final_assembly:
            unoriented = final_assembly.replace('.fasta', '_unoriented.fasta')
            if os.path.exists(unoriented) and stage_fastas:
                stage_fastas = (stage_fastas[:-1] + [(stage_fastas[-1][0], unoriented),
                                                     ('Reoriented (dnaapler)', final_assembly)])

        funnels = ReportMethods.read_funnel(sample, sample_obj, genome_size, cpu)
        assembly_rows = ReportMethods.assembly_funnel(sample, stage_fastas, cpu)
        contigs_long, total_long = ReportMethods.parse_mosdepth_summary(mosdepth_long)
        contigs_short, total_short = ReportMethods.parse_mosdepth_summary(mosdepth_short)

        # Optional QC modules (phase 2): each degrades to a note when its tool or database is missing
        from hashbrown.extra_qc import ExtraQcMethods
        qc_work = os.path.join(report_folder, 'qc', sample)
        Methods.make_folder(qc_work)
        identity, completeness, sanity = {}, {}, {}
        if final_assembly and os.path.exists(final_assembly):
            identity, id_notes = ExtraQcMethods.species_ani(sample, final_assembly, qc_work, cpu)
            notes += id_notes
            completeness, checkm_notes = ExtraQcMethods.checkm2(sample, final_assembly, qc_work, cpu)
            notes += checkm_notes
            busco_result, busco_notes = ExtraQcMethods.busco(sample, final_assembly, qc_work, cpu)
            completeness.update(busco_result)
            notes += busco_notes
            sanity, sanity_notes = ExtraQcMethods.rrna_trna(sample, final_assembly, qc_work, cpu)
            notes += sanity_notes
            qv_map, qv_notes = ExtraQcMethods.merqury_qv(sample, sample_obj, stage_fastas,
                                                         os.path.join(qc_work, 'merqury'), cpu)
            notes += qv_notes
            for row in assembly_rows:
                row['qv'] = qv_map.get(row['stage'], '-')

        # Typing, genes and annotation (phase 3)
        typing, plasmids, genes, annotation = {}, [], [], {}
        if final_assembly and os.path.exists(final_assembly):
            typing, mlst_notes = ExtraQcMethods.mlst(sample, final_assembly)
            notes += mlst_notes
            plasmids, mob_notes = ExtraQcMethods.mob_suite(sample, final_assembly, qc_work, cpu)
            notes += mob_notes
            genes, amr_notes = ExtraQcMethods.amr_virulence(sample, final_assembly, qc_work, cpu)
            notes += amr_notes
            annotation, annotation_notes = ExtraQcMethods.annotate(
                sample, final_assembly, qc_work, cpu, annotate, identity.get('organism', ''))
            notes += annotation_notes
            if annotation.get('gff'):
                links.append((f"Annotation ({annotation['tool']}, GFF3)",
                              os.path.relpath(annotation['gff'], report_folder)))

        # Figures
        read_figs = []
        for funnel in funnels:
            stage_samples = [(row['stage'], *ReportMethods.sample_reads(row['files'][0]))
                             for row in funnel['rows'] if row['files']]
            try:
                svg = ReportMethods.fig_read_distributions(stage_samples)
            except ImportError:
                svg = ''
                notes.append('Figures skipped: matplotlib is not installed.')
            if svg:
                read_figs.append({'title': funnel['title'], 'svg': svg})
        try:
            cov_fig = ReportMethods.fig_contig_coverage(contigs_long, contigs_short)
        except ImportError:
            cov_fig = ''

        # Optional per-stage pycoQC-ng reports
        pycoqc_links = []
        pycoqc_folder = os.path.join(report_folder, 'pycoqc')
        for funnel in funnels:
            kind = 'long' if funnel['title'].startswith('Long') else 'short'
            for row in funnel['rows']:
                tag = f"{kind}_{row['stage'].split(' ')[0]}"
                html = ReportMethods.run_pycoqc(sample, tag, row['files'][0], pycoqc_folder, cpu)
                if html:
                    pycoqc_links.append((f"pycoQC - {funnel['title'].lower()}, {row['stage']}",
                                         os.path.relpath(html, report_folder)))
        if not ReportMethods.pycoqc_available():
            notes.append('pycoQC-ng reports skipped: pycoQC/Fastq_to_seq_summary not on the PATH '
                         '(pip install pycoqc-ng).')

        ReportMethods.write_tables(os.path.join(report_folder, 'tables'), sample, funnels,
                                   assembly_rows, contigs_long, contigs_short,
                                   identity, completeness, sanity,
                                   typing=typing, plasmids=plasmids, genes=genes)

        # Per-contig coverage, long and short side by side
        short_map = {c['contig']: c['mean'] for c in contigs_short}
        cov_rows = [{'contig': c['contig'], 'length': c['length'], 'long': f"{c['mean']:.0f}x",
                     'short': f"{short_map[c['contig']]:.0f}x" if c['contig'] in short_map else '-'}
                    for c in contigs_long]
        cov_rows += [{'contig': c['contig'], 'length': c['length'], 'long': '-',
                      'short': f"{c['mean']:.0f}x"}
                     for c in contigs_short if c['contig'] not in {r['contig'] for r in cov_rows}]
        cov_total = {'length': (total_long or total_short)['length'] if (total_long or total_short) else 0,
                     'long': f"{total_long['mean']:.0f}x" if total_long else '-',
                     'short': f"{total_short['mean']:.0f}x" if total_short else '-'}

        html_file = os.path.join(report_folder, f'{sample}.html')
        with open(html_file, 'w') as f:
            f.write(Template(REPORT_TEMPLATE).render(
                sample=sample, version=__version__, date=date.today().isoformat(),
                sampled=SAMPLED_READS, funnels=funnels, read_figs=read_figs,
                assembly_rows=assembly_rows, cov_rows=cov_rows, cov_total=cov_total,
                cov_fig=cov_fig, identity=identity, completeness=completeness, sanity=sanity,
                typing=typing, plasmids=plasmids, genes=genes, annotation=annotation,
                links=links + pycoqc_links, notes=notes))
        ReportMethods.render_pdf(html_file, notes)

        final = assembly_rows[-1] if assembly_rows else {}
        return {'sample': sample, 'contigs': final.get('contigs', '-'), 'size': final.get('size', '-'),
                'n50': final.get('n50', '-'),
                'long_cov': f"{total_long['mean']:.0f}x" if total_long else '-',
                'short_cov': f"{total_short['mean']:.0f}x" if total_short else '-',
                'st': f"{typing['scheme']} ST{typing['st']}" if typing else '-'}

    @staticmethod
    def build_reports(sample_dict, output_folder, report_folder, mode, cpu, annotate='auto'):
        """One report per sample plus the run-level index.html comparing them."""
        Methods.make_folder(report_folder)
        qc = os.path.join(output_folder, '7_assembly_qc' if mode == 'long' else '3_assembly_qc')
        summaries = []
        for sample, sample_obj in sample_dict.items():
            if mode == 'long':
                stages = [('Assembled (long reads)', '3_assembled_long'),
                          ('Polished - Medaka', '4_polished_long'),
                          ('Polished - short reads', '6_polished_short')]
            else:
                stages = [('Assembled (short reads)', '2_assembled_short')]
            stage_fastas = [(name, os.path.join(output_folder, folder, sample + '.fasta'))
                            for name, folder in stages]
            mosdepth_long = os.path.join(qc, 'mapping_qc_long_reads', sample + '.mosdepth.summary.txt')
            mosdepth_short = os.path.join(qc, 'mapping_qc_short_reads', sample + '.mosdepth.summary.txt')
            links = []
            quast_report = os.path.join(qc, 'quast', 'report.html')
            if os.path.exists(quast_report):
                links.append(('QUAST report', os.path.relpath(quast_report, report_folder)))
            for label, folder in (('Dot plot - raw vs long-polished', 'mummer_long'),
                                  ('Dot plot - long- vs short-polished', 'mummer_short')):
                png = os.path.join(qc, folder, sample + '.png')
                if os.path.exists(png):
                    links.append((label, os.path.relpath(png, report_folder)))
            graph = os.path.join(qc, 'assembly_graphs', sample + '.png')
            if os.path.exists(graph):
                links.append(('Assembly graph (Bandage)', os.path.relpath(graph, report_folder)))
            summaries.append(ReportMethods.build_report(sample, sample_obj, report_folder, stage_fastas,
                                                        mosdepth_long, mosdepth_short, links, cpu,
                                                        annotate))

        with open(os.path.join(report_folder, 'index.html'), 'w') as f:
            f.write(Template(INDEX_TEMPLATE).render(version=__version__, date=date.today().isoformat(),
                                                    summaries=summaries))


STYLE = """
:root { --fg: #24292f; --muted: #57606a; --line: #d0d7de; --bg: #ffffff; --soft: #f6f8fa;
        --accent: #b07219; }
@media (prefers-color-scheme: dark) {
  :root { --fg: #e6edf3; --muted: #8b949e; --line: #30363d; --bg: #0d1117; --soft: #161b22; } }
body { font: 15px/1.5 -apple-system, 'Segoe UI', Roboto, sans-serif; color: var(--fg);
       background: var(--bg); max-width: 960px; margin: 2rem auto; padding: 0 16px; }
h1 { font-size: 1.5rem; border-bottom: 1px solid var(--line); padding-bottom: .4rem; }
h2 { font-size: 1.15rem; margin-top: 2rem; }
table { border-collapse: collapse; margin: .6rem 0; width: 100%; }
th, td { border: 1px solid var(--line); padding: .3rem .6rem; text-align: right; font-size: .88rem; }
th:first-child, td:first-child { text-align: left; }
th { background: var(--soft); }
figure { margin: 1rem 0; overflow-x: auto; }
figure svg { max-width: 100%; height: auto; }
.muted { color: var(--muted); font-size: .85rem; }
.note { background: var(--soft); border-left: 3px solid var(--accent); padding: .4rem .8rem;
        margin: .4rem 0; font-size: .85rem; }
a { color: var(--accent); }
@media print { body { max-width: none; } .note { page-break-inside: avoid; } }
"""

REPORT_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ sample }} - HashBrown report</title><style>""" + STYLE + """</style></head>
<body>
<h1>{{ sample }}</h1>
<p class="muted">HashBrown v{{ version }} &middot; {{ date }}</p>

{% for funnel in funnels %}
<h2>Read QC - {{ funnel.title }}</h2>
<table>
<tr><th>Stage</th><th>Reads</th><th>Bases</th><th>N50</th><th>Mean Q</th><th>%GC</th>
<th>Coverage</th><th>% of raw</th></tr>
{% for r in funnel.rows %}
<tr><td>{{ r.stage }}</td><td>{{ '{:,}'.format(r.reads) }}</td><td>{{ '{:,}'.format(r.bases) }}</td>
<td>{{ '{:,}'.format(r.n50) }}</td><td>{{ r.avg_qual }}</td><td>{{ r.gc }}</td>
<td>{{ r.coverage }}</td><td>{{ r.kept }}</td></tr>
{% endfor %}
</table>
{% endfor %}
{% for fig in read_figs %}
<figure>{{ fig.svg }}<figcaption class="muted">{{ fig.title }}: length and quality by stage
(first {{ sampled }} reads per stage)</figcaption></figure>
{% endfor %}

{% if assembly_rows %}
<h2>Assembly - stage by stage</h2>
<table>
<tr><th>Stage</th><th>Contigs</th><th>Size (bp)</th><th>N50</th><th>%GC</th><th>QV</th></tr>
{% for r in assembly_rows %}
<tr><td>{{ r.stage }}</td><td>{{ '{:,}'.format(r.contigs) }}</td><td>{{ '{:,}'.format(r.size) }}</td>
<td>{{ '{:,}'.format(r.n50) }}</td><td>{{ r.gc }}</td><td>{{ r.qv }}</td></tr>
{% endfor %}
</table>
<p class="muted">QV: Merqury consensus quality from the short-read k-mers (Q30 &asymp; 1 error/kb,
Q50 &asymp; 1 error/100 kb; higher is better).</p>
{% endif %}

{% if identity %}
<h2>Species confirmation</h2>
<table>
<tr><th>Closest genome (mash)</th><th>Mash distance</th><th>Organism</th><th>ANI</th></tr>
<tr><td>{{ identity.hit }}</td><td>{{ identity.distance }}</td>
<td><em>{{ identity.organism }}</em></td><td>{{ identity.ani }}</td></tr>
</table>
<p class="muted">Same species: ANI &ge; 95% against the closest genome.</p>
{% endif %}

{% if completeness or sanity %}
<h2>Completeness &amp; sanity</h2>
<table>
{% if completeness.completeness is defined %}
<tr><td>CheckM2 completeness</td><td>{{ completeness.completeness }}%</td></tr>
<tr><td>CheckM2 contamination</td><td>{{ completeness.contamination }}%</td></tr>
{% endif %}
{% if completeness.summary is defined %}
<tr><td>BUSCO ({{ completeness.lineage }})</td><td>{{ completeness.summary }}</td></tr>
{% endif %}
{% if sanity.rrna is defined %}<tr><td>rRNA genes (barrnap)</td><td>{{ sanity.rrna }}</td></tr>{% endif %}
{% if sanity.trna is defined %}<tr><td>tRNAs (tRNAscan-SE)</td><td>{{ sanity.trna }}</td></tr>{% endif %}
</table>
{% endif %}

{% if typing or plasmids %}
<h2>Typing</h2>
{% if typing %}
<table>
<tr><th>MLST scheme</th><th>ST</th><th>Alleles</th></tr>
<tr><td>{{ typing.scheme }}</td><td>{{ typing.st }}</td><td>{{ typing.alleles }}</td></tr>
</table>
{% endif %}
{% if plasmids %}
<table>
<tr><th>Plasmid (MOB-suite)</th><th>Size (bp)</th><th>Replicon</th><th>Relaxase</th><th>Mobility</th></tr>
{% for p in plasmids %}
<tr><td>{{ p.plasmid }}</td><td>{{ p.size }}</td><td>{{ p.rep }}</td>
<td>{{ p.relaxase }}</td><td>{{ p.mobility }}</td></tr>
{% endfor %}
</table>
{% endif %}
{% endif %}

{% if genes %}
<h2>AMR &amp; virulence genes</h2>
<table>
<tr><th>Gene</th><th>Source</th><th>%Identity</th><th>%Coverage</th><th>Contig</th><th>Class / product</th></tr>
{% for g in genes %}
<tr><td>{{ g.gene }}</td><td>{{ g.source }}</td><td>{{ g.identity }}</td>
<td>{{ g.coverage }}</td><td>{{ g.contig }}</td><td>{{ g.info }}</td></tr>
{% endfor %}
</table>
<p class="muted">AMRFinderPlus: resistance genes and point mutations. VFDB: virulence factors.
PlasmidFinder: plasmid replicons.</p>
{% endif %}

{% if annotation %}
<h2>Annotation</h2>
<p>{{ annotation.tool }}: {{ annotation.summary }}</p>
{% endif %}

{% if cov_rows %}
<h2>Coverage per contig (mosdepth)</h2>
<table>
<tr><th>Contig</th><th>Length (bp)</th><th>Long-read mean</th><th>Short-read mean</th></tr>
{% for c in cov_rows %}
<tr><td>{{ c.contig }}</td><td>{{ '{:,}'.format(c.length) }}</td>
<td>{{ c.long }}</td><td>{{ c.short }}</td></tr>
{% endfor %}
<tr><th>total</th><th>{{ '{:,}'.format(cov_total.length) }}</th>
<th>{{ cov_total.long }}</th><th>{{ cov_total.short }}</th></tr>
</table>
{% if cov_fig %}<figure>{{ cov_fig }}</figure>{% endif %}
{% endif %}

{% if links %}
<h2>Detailed reports</h2>
<ul>{% for label, href in links %}<li><a href="{{ href }}">{{ label }}</a></li>{% endfor %}</ul>
{% endif %}

{% for note in notes %}<div class="note">{{ note }}</div>{% endfor %}
<p class="muted">Tables as TSV: <code>tables/</code></p>
</body></html>
"""

INDEX_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HashBrown run report</title><style>""" + STYLE + """</style></head>
<body>
<h1>HashBrown run</h1>
<p class="muted">HashBrown v{{ version }} &middot; {{ date }}</p>
<table>
<tr><th>Sample</th><th>Contigs</th><th>Size (bp)</th><th>N50</th>
<th>Long-read cov.</th><th>Short-read cov.</th><th>MLST</th></tr>
{% for s in summaries %}
<tr><td><a href="{{ s.sample }}.html">{{ s.sample }}</a></td>
<td>{{ s.contigs }}</td><td>{{ s.size }}</td><td>{{ s.n50 }}</td>
<td>{{ s.long_cov }}</td><td>{{ s.short_cov }}</td><td>{{ s.st }}</td></tr>
{% endfor %}
</table>
</body></html>
"""
