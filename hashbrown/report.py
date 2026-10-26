"""Per-sample HTML/PDF report: the effect of every pipeline stage on the reads and the assembly.

Core tools (seqkit, mosdepth, QUAST) are hard requirements of the pipeline; everything else degrades
gracefully: a missing optional tool (pycoQC-ng, WeasyPrint) is reported as skipped, never an error.
"""

import csv
import gzip
import io
import os
import re
import shutil
import subprocess
from datetime import date, datetime

from jinja2 import Template

from hashbrown import __version__
from hashbrown.methods import Methods

SAMPLED_READS = 10000  # reads per stage used for the distribution figures

# Figure palette: embedded SVGs carry fixed colors, so every one must stay legible on both the light
# (#ffffff) and dark (#0d1117) page. The series order is CVD-checked (deutan/protan/tritan) on both
# surfaces; slot 1 is the HashBrown golden. Ink is the mid gray that reads on both backgrounds.
FIG_SERIES = ('#c87f3a', '#3987e5', '#199e70', '#8a7ee0', '#d55181')
FIG_NEUTRAL = '#898781'  # overflow series / unassigned
FIG_INK = '#767d89'      # axis text, ticks, spines


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
    def read_funnel(sample, sample_obj, genome_size, cpu, mapping_folders=None):
        """One table per read type: raw -> trimmed -> filtered, with the effect of each stage.

        Coverage per stage is the mapped mean depth (that stage's mosdepth summary from the mapping
        QC) when available, and the base-count estimate (bases / genome size, marked with ~) otherwise.
        """
        funnels = []
        mapping_folders = mapping_folders or {}
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

                coverage = f"~{stats['bases'] / genome_size:.0f}x" if genome_size else '-'
                stage_summary = os.path.join(
                    mapping_folders.get(title, ''),
                    f"{sample}_{stage_name.split(' ')[0]}.mosdepth.summary.txt")
                _, total = ReportMethods.parse_mosdepth_summary(stage_summary)
                if total:
                    coverage = f"{total['mean']:.0f}x"

                rows.append({
                    'stage': stage_name,
                    'reads': stats['reads'],
                    'bases': stats['bases'],
                    'n50': stats['n50'],
                    'avg_qual': f"{stats['avg_qual']:.1f}" if stats['avg_qual'] else '-',
                    'gc': f"{stats['gc']:.1f}" if stats['gc'] else '-',
                    'coverage': coverage,
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
    def _thousands(*axes):
        """1,000-style tick labels on the given axes ('x' or 'y' per axis object pair)."""
        from matplotlib.ticker import StrMethodFormatter
        for ax, which in axes:
            getattr(ax, f'{which}axis').set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    @staticmethod
    def _inline_svg(fig):
        """The figure as an SVG string ready to embed in HTML (XML prolog and doctype stripped).

        Saved transparent so the figure sits on the page background in both themes.
        """
        from matplotlib import pyplot as plt
        svg = io.StringIO()
        fig.savefig(svg, format='svg', transparent=True)
        plt.close(fig)
        text = svg.getvalue()
        return text[text.find('<svg'):]

    @staticmethod
    def _style_axes(ax, grid_axis='y'):
        """One look for every figure: recessive spines, ticks and grid in theme-neutral ink."""
        ax.set_facecolor('none')
        ax.spines[['top', 'right']].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set_color(FIG_INK)
            ax.spines[side].set_linewidth(0.8)
        ax.tick_params(colors=FIG_INK, labelcolor=FIG_INK, labelsize=8, width=0.8)
        for axis_label in (ax.xaxis.label, ax.yaxis.label):
            axis_label.set_color(FIG_INK)
            axis_label.set_size(9)
        if grid_axis:
            ax.grid(axis=grid_axis, color=FIG_INK, alpha=0.25, linewidth=0.6)
            ax.set_axisbelow(True)

    @staticmethod
    def _style_legend(ax, **kwargs):
        """Frameless legend with the ink color (figure text keeps one tone in both themes)."""
        legend = ax.legend(frameon=False, fontsize=8, **kwargs)
        for text in legend.get_texts():
            text.set_color(FIG_INK)
        return legend

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
        for i, (name, lengths, quals) in enumerate(stage_samples):
            color = FIG_SERIES[i] if i < len(FIG_SERIES) else FIG_NEUTRAL
            ax_len.hist(lengths, bins=60, histtype='step', linewidth=1.6, label=name, color=color)
            if quals:
                ax_q.hist(quals, bins=40, histtype='step', linewidth=1.6, label=name, color=color)
        ax_len.set_xlabel('Read length (bp)')
        ax_len.set_ylabel('Reads')
        longest = max(max(lengths) for _, lengths, _ in stage_samples)
        shortest = min(min(lengths) for _, lengths, _ in stage_samples)
        if shortest > 0 and longest > 20 * shortest:
            ax_len.set_xscale('log')
        ax_q.set_xlabel('Mean read quality (Phred)')
        ReportMethods._thousands((ax_len, 'x'), (ax_len, 'y'), (ax_q, 'y'))
        for ax in (ax_len, ax_q):
            ReportMethods._style_axes(ax)
        if len(stage_samples) > 1:
            ReportMethods._style_legend(ax_len)
        fig.tight_layout()
        return ReportMethods._inline_svg(fig)

    @staticmethod
    def contig_gc(assembly, cpu):
        """Per-contig length and %GC via seqkit fx2tab; {contig: (length, gc)}."""
        done = subprocess.run(['seqkit', 'fx2tab', '-n', '-i', '-l', '-g', '-j', str(cpu), assembly],
                              capture_output=True, text=True)
        per_contig = {}
        for line in done.stdout.splitlines():
            fields = line.split('\t')
            if len(fields) >= 3:
                try:
                    per_contig[fields[0]] = (int(fields[1]), float(fields[2]))
                except ValueError:
                    continue
        return per_contig

    @staticmethod
    def fig_blob(blob_rows):
        """GC vs coverage per contig, dot area ~ contig length, colored by taxonomy; inline SVG."""
        if not blob_rows:
            return ''
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        organisms = sorted({row['organism'] for row in blob_rows})
        max_length = max(row['length'] for row in blob_rows) or 1
        for i, organism in enumerate(organisms):
            rows = [r for r in blob_rows if r['organism'] == organism]
            color = FIG_SERIES[i] if i < len(FIG_SERIES) else FIG_NEUTRAL
            ax.scatter([r['gc'] for r in rows], [max(r['cov'], 0.1) for r in rows],
                       s=[40 + 360 * r['length'] / max_length for r in rows],
                       color=color, alpha=0.7, edgecolors='none',
                       label=organism[:40] if organism else 'unassigned')
        covs = [max(r['cov'], 0.1) for r in blob_rows]
        if max(covs) > 20 * min(covs):
            ax.set_yscale('log')
        if len(blob_rows) <= 15:
            for r in blob_rows:
                ax.annotate(r['contig'][:15], (r['gc'], max(r['cov'], 0.1)),
                            fontsize=7, color=FIG_INK, xytext=(5, 5), textcoords='offset points')
        ax.set_xlabel('%GC')
        ax.set_ylabel('Mean coverage (x)')
        ReportMethods._thousands((ax, 'y'))
        if any(o for o in organisms):
            ReportMethods._style_legend(ax, loc='best')
        ReportMethods._style_axes(ax, grid_axis='both')
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
                   [long_map.get(n, 0) for n in names], width, label='Long reads', color=FIG_SERIES[0])
        if short_map:
            ax.bar([i + width / 2 for i in x] if long_map else list(x),
                   [short_map.get(n, 0) for n in names], width, label='Short reads', color=FIG_SERIES[1])
        ax.set_xticks(list(x))
        ax.set_xticklabels([n[:20] for n in names], rotation=30, ha='right', fontsize=8)
        ax.set_ylabel('Mean coverage (x)')
        if long_map and short_map:
            ReportMethods._style_legend(ax)
        ReportMethods._style_axes(ax)
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
    def logo_svg():
        """The packaged HashBrown icon as an inline-able SVG string; '' when the asset is missing."""
        try:
            from importlib import resources
            text = (resources.files('hashbrown') / 'logo.svg').read_text()
            return text[text.find('<svg'):]
        except Exception:
            return ''

    @staticmethod
    def _num(value):
        """'99.99%', '100x' or '4,640,168' as a float; None when it is not a number."""
        try:
            return float(str(value).strip().rstrip('%x').replace(',', ''))
        except (TypeError, ValueError):
            return None

    @staticmethod
    def qc_criteria(identity, completeness, sanity, contamination, cov_total):
        """Pass/warn/fail verdicts against the MIMAG standard (Bowers et al. 2017, high-quality draft)
        plus pipeline-sensible read thresholds. Returns (rows, mimag_tier, status_by_key); a metric whose
        module did not run is 'na' (not assessed), never a pass.
        """
        num = ReportMethods._num
        rows = []

        def add(key, metric, measured, threshold, status):
            rows.append({'key': key, 'metric': metric, 'measured': measured,
                         'threshold': threshold, 'status': status})

        comp = num((completeness or {}).get('completeness'))
        add('completeness', 'CheckM2 completeness', f'{comp}%' if comp is not None else '-',
            '> 90% (MIMAG high-quality)',
            'na' if comp is None else 'pass' if comp > 90 else 'warn' if comp >= 50 else 'fail')
        contam = num((completeness or {}).get('contamination'))
        add('checkm2_contamination', 'CheckM2 contamination', f'{contam}%' if contam is not None else '-',
            '< 5% (MIMAG)',
            'na' if contam is None else 'pass' if contam < 5 else 'warn' if contam < 10 else 'fail')

        rrna_counts = None
        if (sanity or {}).get('rrna'):
            rrna_counts = {gene: int(n) for n, gene in re.findall(r'(\d+)x (\w+)', sanity['rrna'])}
        rrna_ok = bool(rrna_counts) and all(rrna_counts.get(g, 0) >= 1 for g in ('5S', '16S', '23S'))
        add('rrna', 'rRNA genes (5S, 16S, 23S)', (sanity or {}).get('rrna', '-'),
            '≥ 1 of each (MIMAG)',
            'na' if rrna_counts is None else 'pass' if rrna_ok else 'fail')

        trna = None
        if (sanity or {}).get('trna'):
            match = re.search(r'(\d+) isotypes', sanity['trna'])
            trna = int(match.group(1)) if match else None
        add('trna', 'tRNA isotypes', (sanity or {}).get('trna', '-'), '≥ 18 (MIMAG)',
            'na' if trna is None else 'pass' if trna >= 18 else 'warn')

        for key, label in (('long', 'Long-read coverage'), ('short', 'Short-read coverage')):
            cov = num((cov_total or {}).get(key))
            if cov is not None:
                add(f'{key}_cov', label, f'{cov:.0f}x', '≥ 30x',
                    'pass' if cov >= 30 else 'warn')

        ani = num((identity or {}).get('ani'))
        add('ani', 'ANI to closest genome', (identity or {}).get('ani', '-'),
            '≥ 95% (same species)',
            'na' if ani is None else 'pass' if ani >= 95 else 'fail')

        verdict = (contamination or {}).get('confindr_status', '')
        snvs = (contamination or {}).get('confindr_snvs', '')
        # The raw SNV count stays visible next to the verdict: ConFindr's contamination
        # threshold differs between genera, so the number lets the reader judge borderline calls
        plural = '' if snvs == '1' else 's'
        measured = f'{verdict} ({snvs} contaminating SNV{plural})' if verdict and snvs != '' else verdict or '-'
        add('confindr', 'ConFindr (intra-species)', measured, 'Clean',
            'na' if not verdict else 'pass' if verdict.lower() == 'clean' else 'fail')

        # The headline: which MIMAG tier the assembly reaches
        if comp is None or contam is None:
            tier = {'label': 'Not assessed', 'status': 'na'}
        elif comp > 90 and contam < 5 and rrna_ok and (trna or 0) >= 18:
            tier = {'label': 'High-quality draft (MIMAG)', 'status': 'pass'}
        elif comp >= 50 and contam < 10:
            tier = {'label': 'Medium-quality draft (MIMAG)', 'status': 'warn'}
        else:
            tier = {'label': 'Low-quality draft (MIMAG)', 'status': 'fail'}
        return rows, tier, {row['key']: row['status'] for row in rows}

    @staticmethod
    def parse_versions(log_file):
        """The 'Tool: version' lines of the run's log.txt as rows for the traceability appendix."""
        rows = []
        if log_file and os.path.exists(log_file):
            with open(log_file) as f:
                for line in f:
                    name, sep, version = line.partition(':')
                    if sep and version.strip():
                        rows.append({'tool': name.strip(), 'version': version.strip()})
        return rows

    @staticmethod
    def write_tables(tables_folder, sample, funnels, assembly_rows, contigs_long, contigs_short,
                     identity=None, completeness=None, sanity=None,
                     typing=None, plasmids=None, genes=None,
                     kmer=None, contamination=None, blob_rows=None):
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
        summary = {**(identity or {}), **(completeness or {}), **(sanity or {}), **(typing or {}),
                   **(kmer or {}), **(contamination or {})}
        if summary:
            write('identity_completeness', list(summary), [summary])
        if plasmids:
            write('plasmids', ['plasmid', 'size', 'rep', 'relaxase', 'mobility'], plasmids)
        if genes:
            write('amr_virulence', ['gene', 'source', 'identity', 'coverage', 'contig', 'info'], genes)
        if blob_rows:
            write('blob', ['contig', 'length', 'gc', 'cov', 'organism'], blob_rows)
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
                     links, cpu, annotate='auto', cmdline='', versions=None):
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

        mapping_folders = {'Long reads': os.path.dirname(mosdepth_long),
                           'Short reads': os.path.dirname(mosdepth_short)}
        funnels = ReportMethods.read_funnel(sample, sample_obj, genome_size, cpu, mapping_folders)
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

        # Contamination and read-based genome profile (phase 4)
        kmer, contamination, blob_rows, blob_fig = {}, {}, [], ''
        if final_assembly and os.path.exists(final_assembly):
            kmer, kmer_notes = ExtraQcMethods.genomescope(
                sample, sample_obj, os.path.join(qc_work, 'merqury', 'reads.meryl'),
                os.path.join(qc_work, 'genomescope'), cpu)
            notes += kmer_notes
            contamination, confindr_notes = ExtraQcMethods.confindr(sample, sample_obj, qc_work, cpu)
            notes += confindr_notes
            cov_map = {c['contig']: c['mean'] for c in (contigs_long or contigs_short)}
            if cov_map:
                taxonomy, tax_notes = ExtraQcMethods.contig_taxonomy(final_assembly, qc_work, cpu)
                notes += tax_notes
                gc_map = ReportMethods.contig_gc(final_assembly, cpu)
                blob_rows = [{'contig': contig, 'length': length, 'gc': gc,
                              'cov': cov_map.get(contig, 0), 'organism': taxonomy.get(contig, '')}
                             for contig, (length, gc) in gc_map.items()]
                try:
                    blob_fig = ReportMethods.fig_blob(blob_rows)
                except ImportError:
                    blob_fig = ''

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
                                   typing=typing, plasmids=plasmids, genes=genes,
                                   kmer=kmer, contamination=contamination, blob_rows=blob_rows)

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

        qc_rows, mimag, qc_status = ReportMethods.qc_criteria(identity, completeness, sanity,
                                                              contamination, cov_total)
        html_file = os.path.join(report_folder, f'{sample}.html')
        with open(html_file, 'w') as f:
            f.write(Template(REPORT_TEMPLATE).render(
                sample=sample, version=__version__, date=date.today().isoformat(),
                timestamp=datetime.now().strftime('%Y-%m-%d %H:%M'), logo=ReportMethods.logo_svg(),
                sampled=SAMPLED_READS, funnels=funnels, read_figs=read_figs,
                assembly_rows=assembly_rows, cov_rows=cov_rows, cov_total=cov_total,
                cov_fig=cov_fig, identity=identity, completeness=completeness, sanity=sanity,
                typing=typing, plasmids=plasmids, genes=genes, annotation=annotation,
                kmer=kmer, contamination=contamination, blob_fig=blob_fig,
                qc_rows=qc_rows, mimag=mimag, qc_status=qc_status,
                cmdline=cmdline, versions=versions or [],
                links=links + pycoqc_links, notes=notes))
        ReportMethods.render_pdf(html_file, notes)

        final = assembly_rows[-1] if assembly_rows else {}
        return {'sample': sample, 'contigs': final.get('contigs', '-'), 'size': final.get('size', '-'),
                'n50': final.get('n50', '-'),
                'long_cov': f"{total_long['mean']:.0f}x" if total_long else '-',
                'short_cov': f"{total_short['mean']:.0f}x" if total_short else '-',
                'st': f"{typing['scheme']} ST{typing['st']}" if typing else '-',
                'mimag': mimag}

    @staticmethod
    def build_reports(sample_dict, output_folder, report_folder, mode, cpu, annotate='auto', cmdline=''):
        """One report per sample plus the run-level index.html comparing them."""
        Methods.make_folder(report_folder)
        qc = os.path.join(output_folder, '7_assembly_qc' if mode == 'long' else '3_assembly_qc')
        versions = ReportMethods.parse_versions(os.path.join(output_folder, 'log.txt'))
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
                                                        annotate, cmdline=cmdline, versions=versions))

        with open(os.path.join(report_folder, 'index.html'), 'w') as f:
            f.write(Template(INDEX_TEMPLATE).render(version=__version__, date=date.today().isoformat(),
                                                    logo=ReportMethods.logo_svg(),
                                                    summaries=summaries))


STYLE = """
:root { --fg: #24292f; --muted: #57606a; --line: #d0d7de; --bg: #ffffff; --soft: #f6f8fa;
        --accent: #b07219; --accent-strong: #8a5a14;
        --good: #1a7f37; --good-soft: #dafbe1; --warn: #9a6700; --warn-soft: #fff3c6;
        --bad: #cf222e; --bad-soft: #ffebe9; }
@media (prefers-color-scheme: dark) {
  :root { --fg: #e6edf3; --muted: #8b949e; --line: #30363d; --bg: #0d1117; --soft: #161b22;
          --accent: #d9983e; --accent-strong: #e3ad62;
          --good: #3fb950; --good-soft: #15261c; --warn: #d29922; --warn-soft: #2b2111;
          --bad: #f85149; --bad-soft: #321a1e; } }
body { font: 15px/1.5 -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
       color: var(--fg); background: var(--bg); max-width: 960px; margin: 0 auto;
       padding: 0 20px 2.5rem; border-top: 4px solid var(--accent); }

/* Masthead: brand row, sample name, verdict and the summary cards */
.masthead { border-bottom: 2px solid var(--line); padding: 1.3rem 0 1.1rem; }
.brand svg { height: 28px; width: auto; vertical-align: -7px; }
.wordmark { font-size: 1.2rem; font-weight: 800; letter-spacing: .02em; margin-left: .4rem; }
.wordmark b { color: var(--accent); }
.ver { font-size: .74rem; font-weight: 600; color: var(--muted); border: 1px solid var(--line);
       border-radius: 999px; padding: .05rem .5rem; margin-left: .5rem; vertical-align: 2px; }
h1 { font-size: 1.65rem; margin: .55rem 0 .1rem; letter-spacing: .01em; }
.verdict { margin: .7rem 0 .2rem; }
.cards { margin: .7rem 0 0; }
.card { display: inline-block; vertical-align: top; border: 1px solid var(--line); border-radius: 8px;
        background: var(--soft); padding: .45rem .75rem .5rem; margin: 0 .4rem .45rem 0;
        min-width: 102px; }
.card .k { font-size: .64rem; font-weight: 700; text-transform: uppercase; letter-spacing: .08em;
           color: var(--muted); }
.card .v { font-size: 1.08rem; font-weight: 700; font-variant-numeric: tabular-nums;
           line-height: 1.35; }
.card .v.sm { font-size: .9rem; line-height: 1.6; }
.card .s { font-size: .73rem; color: var(--muted); }
.st-pass { color: var(--good); } .st-warn { color: var(--warn); } .st-fail { color: var(--bad); }

/* Status chips: always a glyph + word, never color alone */
.chip { display: inline-block; border: 1px solid transparent; border-radius: 999px;
        padding: .02rem .6rem; font-size: .76rem; font-weight: 600; white-space: nowrap; }
.chip.big { font-size: .88rem; padding: .12rem .75rem; }
.status-pass { color: var(--good); background: var(--good-soft); }
.status-warn { color: var(--warn); background: var(--warn-soft); }
.status-fail { color: var(--bad); background: var(--bad-soft); }
.status-na { color: var(--muted); background: var(--soft); border-color: var(--line); }

/* Sections */
h2 { font-size: 1.02rem; margin: 1.9rem 0 .35rem; padding-left: .55rem;
     border-left: 4px solid var(--accent); }
h3 { font-size: .9rem; margin: 1.2rem 0 .2rem; color: var(--muted);
     text-transform: uppercase; letter-spacing: .06em; }
h2 .anchor { color: var(--line); text-decoration: none; margin-left: .4rem; font-weight: 400; }
h2 .anchor:hover { color: var(--accent); }

/* Tables: horizontal rules, zebra stripes, tabular figures */
table { border-collapse: collapse; margin: .6rem 0 1rem; width: auto; max-width: 100%;
        font-variant-numeric: tabular-nums; }
th, td { padding: .32rem 1.1rem .32rem 0; text-align: left; font-size: .86rem; }
th:last-child, td:last-child { padding-right: 0; }
th { font-size: .7rem; font-weight: 700; text-transform: uppercase; letter-spacing: .06em;
     color: var(--muted); border-bottom: 2px solid var(--accent); }
td { border-bottom: 1px solid var(--line); }
tr:nth-child(even) td { background: var(--soft); }
tr.total td { border-top: 2px solid var(--line); border-bottom: none; font-weight: 700;
              background: none; }
td.l, th.l { text-align: left; }  /* kept: harmless now that every column is left-aligned */
td.info { font-size: .78rem; color: var(--muted); }
.versions td { font-size: .8rem; }

/* The AMR table can run to hundreds of rows: scroll it on screen, print it in full */
.scrollbox { max-height: 420px; overflow-y: auto; border: 1px solid var(--line);
             border-radius: 6px; }
.scrollbox table { margin: 0; }
.scrollbox tr:first-child th { position: sticky; top: 0; background: var(--bg); }

figure { margin: 1.1rem 0; overflow-x: auto; }
figure svg { max-width: 100%; height: auto; }
.muted { color: var(--muted); font-size: .85rem; }
.note { background: var(--soft); border-left: 3px solid var(--accent); padding: .4rem .8rem;
        margin: .4rem 0; font-size: .85rem; }
a { color: var(--accent); }
code { background: var(--soft); border-radius: 4px; padding: .05rem .3rem; font-size: .85em; }
.cmd { font-family: ui-monospace, 'SFMono-Regular', Consolas, 'Liberation Mono', monospace;
       word-break: break-all; }

@media print {
  body { max-width: none; border-top: none; padding: 0; }
  .scrollbox { max-height: none; overflow: visible; border: none; border-radius: 0; }
  .scrollbox tr:first-child th { position: static; }
  .scrollbox table { margin: .6rem 0 1rem; }
  tr, .note, .card, figure { page-break-inside: avoid; }
  h1, h2, h3 { page-break-after: avoid; }
  h2 .anchor { display: none; }
}
"""

# Shared template pieces: the brand row and the status chip (glyph + word, never color alone)
MACROS = """{% macro chip(status, label=none, big=false) -%}
<span class="chip{{ ' big' if big }} status-{{ status }}">{{
  {'pass': '✓', 'warn': '!', 'fail': '✗', 'na': '–'}[status] }}&nbsp;{{
  label if label is not none else
  {'pass': 'Pass', 'warn': 'Warn', 'fail': 'Fail', 'na': 'Not assessed'}[status] }}</span>
{%- endmacro %}
{% macro brand(version, logo) -%}
<div class="brand">{{ logo }}<span class="wordmark">Hash<b>Brown</b></span>
<span class="ver">v{{ version }}</span></div>
{%- endmacro %}"""

REPORT_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ sample }} - HashBrown report</title><style>""" + STYLE + """</style>
<style>
@page { margin: 16mm 14mm 18mm;
  @bottom-left { content: "{{ sample }} · HashBrown v{{ version }} · {{ date }}";
                 font-size: 8pt; color: #8b949e; }
  @bottom-right { content: "Page " counter(page) " of " counter(pages);
                  font-size: 8pt; color: #8b949e; } }
</style></head>
<body>""" + MACROS + """
<header class="masthead">
{{ brand(version, logo) }}
<h1>{{ sample }}</h1>
<p class="muted">Assembly report &middot; {{ date }}</p>
{% if mimag %}<p class="verdict">{{ chip(mimag.status, mimag.label, big=true) }}</p>{% endif %}
<div class="cards">
{% if assembly_rows %}{% set final = assembly_rows[-1] %}
<div class="card"><div class="k">Genome size</div>
<div class="v">{{ '{:,}'.format(final.size) }}</div><div class="s">bp</div></div>
<div class="card"><div class="k">Contigs</div>
<div class="v">{{ '{:,}'.format(final.contigs) }}</div><div class="s">%GC {{ final.gc }}</div></div>
<div class="card"><div class="k">N50</div>
<div class="v">{{ '{:,}'.format(final.n50) }}</div><div class="s">bp</div></div>
{% endif %}
{% if cov_total and (cov_total.long != '-' or cov_total.short != '-') %}
<div class="card"><div class="k">Coverage</div>
<div class="v">{{ cov_total.long }} / {{ cov_total.short }}</div>
<div class="s">long / short reads</div></div>
{% endif %}
{% if identity.organism is defined and identity.organism != '-' %}
<div class="card"><div class="k">Species</div>
<div class="v sm"><em>{{ identity.organism|truncate(34) }}</em></div>
<div class="s st-{{ qc_status.get('ani', 'na') }}">ANI {{ identity.ani }}</div></div>
{% endif %}
{% if typing %}
<div class="card"><div class="k">MLST</div>
<div class="v">ST{{ typing.st }}</div><div class="s">{{ typing.scheme }}</div></div>
{% endif %}
{% if completeness.completeness is defined %}
<div class="card"><div class="k">CheckM2</div>
<div class="v st-{{ qc_status.get('completeness', 'na') }}">{{ completeness.completeness }}%</div>
<div class="s st-{{ qc_status.get('checkm2_contamination', 'na') }}">
contam. {{ completeness.contamination }}%</div></div>
{% endif %}
{% if contamination.confindr_status is defined %}
<div class="card"><div class="k">ConFindr</div>
<div class="v sm st-{{ qc_status.get('confindr', 'na') }}">{{ contamination.confindr_status }}{%
 if contamination.confindr_snvs is defined %} ({{ contamination.confindr_snvs }} SNV{{
 '' if contamination.confindr_snvs == '1' else 's' }}){% endif %}</div>
<div class="s">{{ contamination.confindr_genus }}</div></div>
{% endif %}
</div>
</header>

{% if qc_rows %}
<h2 id="qc-criteria">QC criteria<a class="anchor" href="#qc-criteria">#</a></h2>
<table>
<thead>
<tr><th>Metric</th><th>Measured</th><th class="l">Threshold</th><th>Status</th></tr>
</thead>
<tbody>
{% for r in qc_rows %}
<tr><td>{{ r.metric }}</td><td>{{ r.measured }}</td><td class="l">{{ r.threshold }}</td>
<td>{{ chip(r.status) }}</td></tr>
{% endfor %}
</tbody>
</table>
<p class="muted">Thresholds follow the MIMAG high-quality draft standard (Bowers et al. 2017,
doi:10.1038/nbt.3893) plus pipeline read-coverage minimums. A metric whose module did not run is
reported as not assessed, never as a pass.</p>
{% endif %}

{% for funnel in funnels %}
<h2 id="read-qc-{{ funnel.title|lower|replace(' ', '-') }}">Read QC - {{ funnel.title }}<a
 class="anchor" href="#read-qc-{{ funnel.title|lower|replace(' ', '-') }}">#</a></h2>
<table>
<thead>
<tr><th>Stage</th><th>Reads</th><th>Bases</th><th>N50</th><th>Mean Q</th><th>%GC</th>
<th>Coverage</th><th>% of raw</th></tr>
</thead>
<tbody>
{% for r in funnel.rows %}
<tr><td>{{ r.stage }}</td><td>{{ '{:,}'.format(r.reads) }}</td><td>{{ '{:,}'.format(r.bases) }}</td>
<td>{{ '{:,}'.format(r.n50) }}</td><td>{{ r.avg_qual }}</td><td>{{ r.gc }}</td>
<td>{{ r.coverage }}</td><td>{{ r.kept }}</td></tr>
{% endfor %}
</tbody>
</table>
<p class="muted">Coverage is the mapped mean depth of that stage's reads on the final assembly
(mosdepth); a ~ marks a base-count estimate where the stage was not mapped.</p>
{% endfor %}
{% for fig in read_figs %}
<figure>{{ fig.svg }}<figcaption class="muted">{{ fig.title }}: length and quality by stage
(first {{ sampled }} reads per stage)</figcaption></figure>
{% endfor %}

{% if assembly_rows %}
<h2 id="assembly">Assembly - stage by stage<a class="anchor" href="#assembly">#</a></h2>
<table>
<thead>
<tr><th>Stage</th><th>Contigs</th><th>Size (bp)</th><th>N50</th><th>%GC</th><th>QV</th></tr>
</thead>
<tbody>
{% for r in assembly_rows %}
<tr><td>{{ r.stage }}</td><td>{{ '{:,}'.format(r.contigs) }}</td><td>{{ '{:,}'.format(r.size) }}</td>
<td>{{ '{:,}'.format(r.n50) }}</td><td>{{ r.gc }}</td><td>{{ r.qv }}</td></tr>
{% endfor %}
</tbody>
</table>
<p class="muted">QV: Merqury consensus quality from the short-read k-mers (Q30 &asymp; 1 error/kb,
Q50 &asymp; 1 error/100 kb; higher is better).</p>
{% endif %}

{% if identity %}
<h2 id="species">Species confirmation<a class="anchor" href="#species">#</a></h2>
<table>
<thead>
<tr><th>Closest genome (mash)</th><th>Mash distance</th><th>Organism</th><th>ANI</th></tr>
</thead>
<tbody>
<tr><td>{{ identity.hit }}</td><td>{{ identity.distance }}</td>
<td class="l"><em>{{ identity.organism }}</em></td><td>{{ identity.ani }}</td></tr>
</tbody>
</table>
<p class="muted">Same species: ANI &ge; 95% against the closest genome.</p>
{% endif %}

{% if completeness or sanity %}
<h2 id="completeness">Completeness &amp; sanity<a class="anchor" href="#completeness">#</a></h2>
<table>
<thead>
<tr><th>Metric</th><th>Value</th></tr>
</thead>
<tbody>
{% if completeness.completeness is defined %}
<tr><td>CheckM2 completeness</td><td>{{ completeness.completeness }}%</td></tr>
<tr><td>CheckM2 contamination</td><td>{{ completeness.contamination }}%</td></tr>
{% endif %}
{% if completeness.summary is defined %}
<tr><td>BUSCO ({{ completeness.lineage }})</td><td>{{ completeness.summary }}</td></tr>
{% endif %}
{% if sanity.rrna is defined %}<tr><td>rRNA genes (barrnap)</td><td>{{ sanity.rrna }}</td></tr>{% endif %}
{% if sanity.trna is defined %}<tr><td>tRNAs (tRNAscan-SE)</td><td>{{ sanity.trna }}</td></tr>{% endif %}
</tbody>
</table>
{% endif %}

{% if typing or plasmids %}
<h2 id="typing">Typing<a class="anchor" href="#typing">#</a></h2>
{% if typing %}
<table>
<thead>
<tr><th>MLST scheme</th><th>ST</th><th>Alleles</th></tr>
</thead>
<tbody>
<tr><td>{{ typing.scheme }}</td><td>{{ typing.st }}</td><td>{{ typing.alleles }}</td></tr>
</tbody>
</table>
{% endif %}
{% if plasmids %}
<table>
<thead>
<tr><th>Plasmid (MOB-suite)</th><th>Size (bp)</th><th>Replicon</th><th>Relaxase</th><th>Mobility</th></tr>
</thead>
<tbody>
{% for p in plasmids %}
<tr><td>{{ p.plasmid }}</td><td>{{ p.size }}</td><td>{{ p.rep }}</td>
<td>{{ p.relaxase }}</td><td>{{ p.mobility }}</td></tr>
{% endfor %}
</tbody>
</table>
{% endif %}
{% endif %}

{% if genes %}
<h2 id="amr">AMR &amp; virulence genes <span class="muted">({{ genes|length }})</span><a
 class="anchor" href="#amr">#</a></h2>
<div class="scrollbox">
<table>
<thead>
<tr><th>Gene</th><th>Source</th><th>%Identity</th><th>%Coverage</th><th>Contig</th><th>Class / product</th></tr>
</thead>
<tbody>
{% for g in genes %}
<tr><td>{{ g.gene }}</td><td>{{ g.source }}</td><td>{{ g.identity }}</td>
<td>{{ g.coverage }}</td><td>{{ g.contig }}</td><td class="info">{{ g.info }}</td></tr>
{% endfor %}
</tbody>
</table>
</div>
<p class="muted">AMRFinderPlus: resistance genes and point mutations. VFDB: virulence factors.
PlasmidFinder: plasmid replicons.</p>
{% endif %}

{% if annotation %}
<h2 id="annotation">Annotation<a class="anchor" href="#annotation">#</a></h2>
<p>{{ annotation.tool }}: {{ annotation.summary }}</p>
{% endif %}

{% if kmer %}
<h2 id="genomescope">Read-based genome profile (GenomeScope2)<a class="anchor"
 href="#genomescope">#</a></h2>
<table>
<thead>
<tr><th>Property</th><th>Value</th></tr>
</thead>
<tbody>
{% if kmer.genome_size_estimate is defined %}
<tr><td>Estimated haploid genome size</td><td>{{ kmer.genome_size_estimate }}</td></tr>
{% endif %}
{% if kmer.heterozygosity is defined %}
<tr><td>Heterozygosity (should be ~0 for a pure isolate)</td><td>{{ kmer.heterozygosity }}</td></tr>
{% endif %}
{% if kmer.read_error_rate is defined %}
<tr><td>Read error rate</td><td>{{ kmer.read_error_rate }}</td></tr>
{% endif %}
</tbody>
</table>
<p class="muted">From the short-read 21-mers, independent of the assembly: a size estimate far from the
assembly size, or non-zero heterozygosity, points at contamination or a mixed culture.</p>
{% endif %}

{% if contamination or blob_fig %}
<h2 id="contamination">Contamination<a class="anchor" href="#contamination">#</a></h2>
{% if contamination %}
<table>
<thead>
<tr><th>ConFindr verdict</th><th>Genus</th><th>Contaminating SNVs</th></tr>
</thead>
<tbody>
<tr><td>{{ contamination.confindr_status }}</td><td>{{ contamination.confindr_genus }}</td>
<td>{{ contamination.confindr_snvs }}</td></tr>
</tbody>
</table>
{% endif %}
{% if blob_fig %}
<figure>{{ blob_fig }}<figcaption class="muted">Blob plot: every contig by %GC and coverage (dot area ~
length, color = closest RefSeq genome). Contigs far from the main blob, or colored as another species,
are contamination candidates.</figcaption></figure>
{% endif %}
{% endif %}

{% if cov_rows %}
<h2 id="coverage">Coverage per contig (mosdepth)<a class="anchor" href="#coverage">#</a></h2>
<table>
<thead>
<tr><th>Contig</th><th>Length (bp)</th><th>Long-read mean</th><th>Short-read mean</th></tr>
</thead>
<tbody>
{% for c in cov_rows %}
<tr><td>{{ c.contig }}</td><td>{{ '{:,}'.format(c.length) }}</td>
<td>{{ c.long }}</td><td>{{ c.short }}</td></tr>
{% endfor %}
<tr class="total"><td>total</td><td>{{ '{:,}'.format(cov_total.length) }}</td>
<td>{{ cov_total.long }}</td><td>{{ cov_total.short }}</td></tr>
</tbody>
</table>
{% if cov_fig %}<figure>{{ cov_fig }}</figure>{% endif %}
{% endif %}

{% if links %}
<h2 id="links">Detailed reports<a class="anchor" href="#links">#</a></h2>
<ul>{% for label, href in links %}<li><a href="{{ href }}">{{ label }}</a></li>{% endfor %}</ul>
{% endif %}

{% if versions or cmdline %}
<h2 id="traceability">Traceability<a class="anchor" href="#traceability">#</a></h2>
<table>
<thead>
<tr><th>Item</th><th>Value</th></tr>
</thead>
<tbody>
<tr><td>HashBrown</td><td class="l">v{{ version }}</td></tr>
<tr><td>Report generated</td><td class="l">{{ timestamp }}</td></tr>
{% if cmdline %}
<tr><td>Command line</td><td class="l"><code class="cmd">{{ cmdline }}</code></td></tr>
{% endif %}
</tbody>
</table>
{% if versions %}
<h3>Software versions</h3>
<table class="versions">
<thead>
<tr><th>Software</th><th>Version</th></tr>
</thead>
<tbody>
{% for v in versions %}<tr><td>{{ v.tool }}</td><td class="l">{{ v.version }}</td></tr>{% endfor %}
</tbody>
</table>
{% endif %}
{% endif %}

{% for note in notes %}<div class="note">{{ note }}</div>{% endfor %}
<p class="muted">Tables as TSV: <code>tables/</code></p>
</body></html>
"""

INDEX_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HashBrown run report</title><style>""" + STYLE + """</style></head>
<body>""" + MACROS + """
<header class="masthead">
{{ brand(version, logo) }}
<h1>HashBrown run</h1>
<p class="muted">{{ summaries|length }} sample{{ 's' if summaries|length != 1 }} &middot;
{{ date }}</p>
</header>
<table>
<thead>
<tr><th>Sample</th><th>Contigs</th><th>Size (bp)</th><th>N50</th>
<th>Long-read cov.</th><th>Short-read cov.</th><th>MLST</th><th>Verdict</th></tr>
</thead>
<tbody>
{% for s in summaries %}
<tr><td><a href="{{ s.sample }}.html">{{ s.sample }}</a></td>
<td>{{ s.contigs }}</td><td>{{ s.size }}</td><td>{{ s.n50 }}</td>
<td>{{ s.long_cov }}</td><td>{{ s.short_cov }}</td><td>{{ s.st }}</td>
<td class="l">{% if s.mimag %}{{ chip(s.mimag.status, s.mimag.label) }}{% endif %}</td></tr>
{% endfor %}
</tbody>
</table>
</body></html>
"""
