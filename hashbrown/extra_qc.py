"""Optional assembly QC modules (phase 2 of docs/reporting_plan.md).

Every module runs only when its tool (and database, where one is needed) is available; otherwise it
returns empty results plus a note saying how to enable it. Nothing here may ever fail the pipeline.
"""

import contextlib
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import zipfile

from hashbrown.methods import Methods


def db_dir():
    """Database folder for the optional modules; $HASHBROWN_DB overrides the default."""
    return os.environ.get('HASHBROWN_DB', os.path.join(os.path.expanduser('~'), '.hashbrown', 'db'))


def _run(cmd, **kwargs):
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kwargs)


def _missing(*tools):
    return [t for t in tools if not shutil.which(t)]


class ExtraQcMethods:
    # ------------------------------------------------------------------ dnaapler
    @staticmethod
    def reorient_dnaapler(sample, sample_obj, cpu):
        """Reorient the final assembly in place (chromosome to dnaA, plasmids to rep) with dnaapler.

        The pre-reorientation sequence is kept next to it as <sample>_unoriented.fasta.
        Returns True when the assembly was reoriented.
        """
        final = sample_obj.best_assembly()
        if not final or not os.path.exists(final):
            return False
        work = os.path.join(os.path.dirname(final), sample + '_dnaapler')
        done = _run(['dnaapler', 'all', '-i', final, '-o', work, '-p', sample,
                     '-t', str(cpu), '-f'])
        reoriented = glob.glob(os.path.join(work, '*reoriented.fasta'))
        if done.returncode == 0 and reoriented and os.path.getsize(reoriented[0]):
            unoriented = os.path.join(os.path.dirname(final),
                                      os.path.basename(final).replace('.fasta', '_unoriented.fasta'))
            shutil.move(final, unoriented)
            Methods.format_fasta(reoriented[0], final)
            sample_obj.assembly.unoriented = unoriented
            shutil.rmtree(work, ignore_errors=True)
            return True
        shutil.rmtree(work, ignore_errors=True)
        return False

    @staticmethod
    def reorient_parallel(sample_dict, done_flag, cpu):
        """Reorient every sample's final assembly; resume-safe via the done flag."""
        if os.path.exists(done_flag):
            print('\tSkipping reorienting. Already done.')
            return
        if _missing('dnaapler'):
            print('\tSkipped: dnaapler is not installed (conda install dnaapler).')
            return
        for sample, sample_obj in sample_dict.items():
            if ExtraQcMethods.reorient_dnaapler(sample, sample_obj, cpu):
                print(f'\t{sample}')
        open(done_flag, 'w').close()

    # ------------------------------------------------------------------ CheckM2
    @staticmethod
    def checkm2(sample, assembly, work, cpu):
        """Completeness and contamination; returns ({'completeness':, 'contamination':}, notes)."""
        if _missing('checkm2'):
            return {}, ['CheckM2 skipped: not installed (conda install -c bioconda checkm2, '
                        'then hashbrown db checkm2).']
        out = os.path.join(work, 'checkm2')
        report = os.path.join(out, 'quality_report.tsv')
        if not os.path.exists(report):
            db = glob.glob(os.path.join(db_dir(), 'checkm2', '**', '*.dmnd'), recursive=True)
            cmd = ['checkm2', 'predict', '-i', assembly, '-o', out, '-t', str(cpu), '--force']
            if db:
                cmd += ['--database_path', db[0]]
            _run(cmd)
        if not os.path.exists(report):
            return {}, ['CheckM2 failed: is its database installed? (hashbrown db checkm2)']
        with open(report) as f:
            header = f.readline().rstrip('\n').split('\t')
            fields = f.readline().rstrip('\n').split('\t')
        row = dict(zip(header, fields, strict=False))
        return {'completeness': row.get('Completeness', '-'),
                'contamination': row.get('Contamination', '-')}, []

    # ------------------------------------------------------------------ BUSCO
    @staticmethod
    def busco(sample, assembly, work, cpu):
        """BUSCO with prokaryote auto-lineage; returns ({'summary':, 'lineage':}, notes)."""
        if _missing('busco'):
            return {}, ['BUSCO skipped: not installed (conda install -c bioconda busco).']
        out_name = f'{sample}_busco'
        summaries = glob.glob(os.path.join(work, out_name, 'short_summary.*.txt'))
        if not summaries:
            _run(['busco', '-i', assembly, '-m', 'genome', '--auto-lineage-prok',
                  '-o', out_name, '--out_path', work, '-c', str(cpu), '-f', '-q'])
            summaries = glob.glob(os.path.join(work, out_name, 'short_summary.*.txt'))
        if not summaries:
            return {}, ['BUSCO failed: see its logs (it downloads lineages on first use; '
                        'an offline machine needs busco --download prokaryota first).']
        summary_file = sorted(summaries, key=lambda p: 'specific' not in p)[0]
        summary, lineage = '-', '-'
        with open(summary_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('C:'):
                    summary = line
                elif 'dataset is:' in line:
                    lineage = line.split('dataset is:')[1].split('(')[0].strip()
        return {'summary': summary, 'lineage': lineage}, []

    # ------------------------------------------------------------------ Merqury
    @staticmethod
    def merqury_qv(sample, sample_obj, stage_fastas, work, cpu):
        """Reference-free consensus QV per assembly stage from short-read k-mers (meryl + merqury).

        Returns ({stage_name: qv}, notes).
        """
        reads = sample_obj.best_illumina()
        if not reads:
            return {}, []
        if _missing('meryl', 'merqury.sh'):
            return {}, ['Merqury QV skipped: not installed (conda install -c bioconda merqury).']
        meryl_db = os.path.join(work, 'reads.meryl')
        if not ExtraQcMethods._build_meryl_db(reads, meryl_db, cpu):
            return {}, ['Merqury QV failed: meryl could not count the short-read k-mers.']

        qv_per_stage = {}
        for stage_name, fasta in stage_fastas:
            if not fasta or not os.path.exists(fasta):
                continue
            tag = re.sub(r'\W+', '_', stage_name.lower()).strip('_')
            stage_dir = os.path.join(work, tag)
            qv_file = os.path.join(stage_dir, f'{tag}.qv')
            if not os.path.exists(qv_file):
                Methods.make_folder(stage_dir)
                _run(['merqury.sh', os.path.abspath(meryl_db), os.path.abspath(fasta), tag],
                     cwd=stage_dir)
            if os.path.exists(qv_file):
                with open(qv_file) as f:
                    fields = f.readline().split('\t')
                if len(fields) >= 4:
                    qv_per_stage[stage_name] = fields[3]
        return qv_per_stage, [] if qv_per_stage else ['Merqury QV failed: no .qv output produced.']

    @staticmethod
    def _build_meryl_db(reads, meryl_db, cpu):
        """21-mer counts of the short reads, shared by Merqury and GenomeScope2; True when usable."""
        if os.path.isdir(meryl_db):
            return True
        Methods.make_folder(os.path.dirname(meryl_db))
        parts = []
        for i, fastq in enumerate(reads):
            part = os.path.join(os.path.dirname(meryl_db), f'reads_{i}.meryl')
            _run(['meryl', 'count', 'k=21', f'threads={cpu}', fastq, 'output', part])
            parts.append(part)
        _run(['meryl', 'union-sum', 'output', meryl_db] + parts)
        for part in parts:
            shutil.rmtree(part, ignore_errors=True)
        return os.path.isdir(meryl_db)

    # ------------------------------------------------------------------ GenomeScope2
    @staticmethod
    def genomescope(sample, sample_obj, meryl_db, work, cpu):
        """Read-based genome profile from the short-read 21-mers; returns (result dict, notes)."""
        reads = sample_obj.best_illumina()
        if not reads:
            return {}, []
        if _missing('genomescope2'):
            return {}, ['GenomeScope2 skipped: not installed (conda install genomescope2).']
        if _missing('meryl'):
            return {}, ['GenomeScope2 skipped: it needs meryl for the k-mer histogram '
                        '(conda install merqury).']
        Methods.make_folder(work)
        if not ExtraQcMethods._build_meryl_db(reads, meryl_db, cpu):
            return {}, ['GenomeScope2 failed: meryl could not count the short-read k-mers.']
        histogram = os.path.join(work, 'reads.histo')
        if not os.path.exists(histogram):
            done = subprocess.run(['meryl', 'histogram', meryl_db], capture_output=True, text=True)
            if done.returncode != 0 or not done.stdout:
                return {}, ['GenomeScope2 failed: meryl histogram produced nothing.']
            with open(histogram, 'w') as f:
                f.write(done.stdout)
        out = os.path.join(work, 'genomescope')
        summary_file = os.path.join(out, 'summary.txt')
        if not os.path.exists(summary_file):
            _run(['genomescope2', '-i', histogram, '-o', out, '-k', '21', '-p', '1'])
        if not os.path.exists(summary_file):
            return {}, ['GenomeScope2 failed: see its output in the report qc folder.']
        results = {}
        wanted = {'Genome Haploid Length': 'genome_size_estimate',
                  'Heterozygous (ab)': 'heterozygosity',
                  'Read Error Rate': 'read_error_rate'}
        with open(summary_file) as f:
            for line in f:
                for label, key in wanted.items():
                    if line.startswith(label):
                        values = re.split(r'\s{2,}', line.strip())
                        if len(values) >= 2:
                            results[key] = ' - '.join(values[1:3]) if len(values) > 2 else values[1]
        return results, []

    # ------------------------------------------------------------------ ConFindr
    @staticmethod
    def confindr(sample, sample_obj, work, cpu):
        """Intra-species contamination from the short reads; returns (result dict, notes)."""
        reads = sample_obj.best_illumina()
        if not reads:
            return {}, []
        exe = shutil.which('confindr.py') or shutil.which('confindr')
        if not exe:
            return {}, ['ConFindr skipped: not installed (conda install confindr).']
        out = os.path.join(work, 'confindr')
        report = os.path.join(out, 'confindr_report.csv')
        if not os.path.exists(report):
            reads_dir = os.path.join(work, 'confindr_reads')
            Methods.make_folder(reads_dir)
            for fastq in reads:
                link = os.path.join(reads_dir, os.path.basename(fastq))
                if not os.path.exists(link):
                    os.symlink(os.path.realpath(fastq), link)
            cmd = [exe, '-i', reads_dir, '-o', out, '-t', str(cpu)]
            database = os.path.join(db_dir(), 'confindr')
            if os.path.isdir(database):
                cmd += ['-d', database]
            _run(cmd)
        if not os.path.exists(report):
            return {}, ['ConFindr failed: its rMLST database may be missing (hashbrown db confindr; '
                        'full coverage needs a free PubMLST account, see the ConFindr docs).']
        with open(report) as f:
            rows = list(csv.DictReader(f))
        if not rows:
            return {}, ['ConFindr produced an empty report.']
        contaminated = any(r.get('ContamStatus', '').strip().lower() == 'true' for r in rows)
        return {'confindr_status': 'Contaminated' if contaminated else 'Clean',
                'confindr_genus': ';'.join(sorted({r.get('Genus') or '-' for r in rows})),
                'confindr_snvs': ';'.join(r.get('NumContamSNVs') or '-' for r in rows)}, []

    # ------------------------------------------------------------------ per-contig taxonomy
    @staticmethod
    def contig_taxonomy(assembly, work, cpu):
        """Mash top hit per contig (RefSeq sketch), named via NCBI datasets when online.

        Returns ({contig: organism or accession}, notes). Empty when mash or its database is absent —
        the blob plot then draws uncolored.
        """
        if _missing('mash'):
            return {}, []
        sketches = glob.glob(os.path.join(db_dir(), 'mash', '*.msh'))
        if not sketches:
            return {}, []
        done = subprocess.run(['mash', 'dist', '-p', str(cpu), '-i', sketches[0], assembly],
                              capture_output=True, text=True)
        best = {}
        for line in done.stdout.splitlines():
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            contig, distance = fields[1], float(fields[2])
            if contig not in best or distance < best[contig][1]:
                best[contig] = (fields[0], distance)
        if not best:
            return {}, []

        # Name the few unique accessions (cached, so a rebuilt report is offline-safe)
        cache_file = os.path.join(work, 'organisms.json')
        organisms = {}
        if os.path.exists(cache_file):
            with open(cache_file) as f:
                organisms = json.load(f)
        accessions = {}
        for contig, (hit, _) in best.items():
            match = re.search(r'(GC[AF]_\d+\.\d+)', hit)
            accessions[contig] = match.group(1) if match else os.path.basename(hit)
        if shutil.which('datasets'):
            for accession in set(accessions.values()) - set(organisms):
                if not accession.startswith(('GCA_', 'GCF_')):
                    continue
                summary = subprocess.run(['datasets', 'summary', 'genome', 'accession', accession],
                                         capture_output=True, text=True)
                with contextlib.suppress(json.JSONDecodeError, KeyError, IndexError):
                    organisms[accession] = \
                        json.loads(summary.stdout)['reports'][0]['organism']['organism_name']
            with open(cache_file, 'w') as f:
                json.dump(organisms, f)
        return {contig: organisms.get(acc, acc) for contig, acc in accessions.items()}, []

    # ------------------------------------------------------------------ rRNA / tRNA
    @staticmethod
    def rrna_trna(sample, assembly, work, cpu):
        """Fast annotation sanity check: rRNA genes (barrnap) and tRNA isotypes (tRNAscan-SE)."""
        results, notes = {}, []
        if _missing('barrnap'):
            notes.append('rRNA counts skipped: barrnap is not installed (conda install barrnap).')
        else:
            done = subprocess.run(['barrnap', '--threads', str(cpu), assembly],
                                  capture_output=True, text=True)
            counts = {}
            for line in done.stdout.splitlines():
                match = re.search(r'Name=(\w+?)_rRNA', line)
                if match:
                    counts[match.group(1)] = counts.get(match.group(1), 0) + 1
            if done.returncode == 0:
                results['rrna'] = ', '.join(f'{counts.get(gene, 0)}x {gene}'
                                            for gene in ('5S', '16S', '23S'))
        if _missing('tRNAscan-SE'):
            notes.append('tRNA counts skipped: tRNAscan-SE is not installed '
                         '(conda install trnascan-se).')
        else:
            trna_out = os.path.join(work, 'trnascan.txt')
            if not os.path.exists(trna_out):
                _run(['tRNAscan-SE', '-B', '-q', '-o', trna_out, assembly])
            if os.path.exists(trna_out):
                with open(trna_out) as f:
                    rows = [line.split('\t') for line in f.readlines()[3:] if '\t' in line]
                isotypes = {fields[4].strip() for fields in rows if len(fields) > 4}
                results['trna'] = f'{len(rows)} tRNAs, {len(isotypes)} isotypes'
        return results, notes

    # ------------------------------------------------------------------ species ID + ANI
    @staticmethod
    def species_ani(sample, assembly, work, cpu):
        """Mash species ID against the local sketch database, then ANI against the top hit.

        Returns ({'hit':, 'distance':, 'organism':, 'ani':}, notes).
        """
        notes = []
        if _missing('mash'):
            return {}, ['Species ID skipped: mash is not installed (conda install mash).']
        sketches = glob.glob(os.path.join(db_dir(), 'mash', '*.msh'))
        if not sketches:
            return {}, ['Species ID skipped: no mash sketch database (run hashbrown db mash).']
        done = subprocess.run(['mash', 'dist', '-p', str(cpu), sketches[0], assembly],
                              capture_output=True, text=True)
        hits = [line.split('\t') for line in done.stdout.splitlines() if '\t' in line]
        if done.returncode != 0 or not hits:
            return {}, ['Species ID failed: mash dist produced no hits.']
        top = min(hits, key=lambda fields: float(fields[2]))
        results = {'hit': os.path.basename(top[0]), 'distance': top[2], 'organism': '-', 'ani': '-'}

        match = re.search(r'(GC[AF]_\d+\.\d+)', top[0])
        if not match:
            return results, notes + ['ANI skipped: no assembly accession in the mash hit name.']
        accession = match.group(1)
        if _missing('datasets'):
            return results, notes + ['ANI skipped: NCBI datasets is not installed '
                                     '(conda install ncbi-datasets-cli).']

        summary = subprocess.run(['datasets', 'summary', 'genome', 'accession', accession],
                                 capture_output=True, text=True)
        try:
            report = json.loads(summary.stdout)['reports'][0]
            results['organism'] = report['organism']['organism_name']
        except (json.JSONDecodeError, KeyError, IndexError):
            notes.append(f'Could not fetch the organism name of {accession} (offline?).')

        ref_fasta = ExtraQcMethods._fetch_genome(accession, work)
        if not ref_fasta:
            return results, notes + [f'ANI skipped: could not download {accession} (offline?).']
        ani = ExtraQcMethods._ani(assembly, ref_fasta, work, cpu)
        if ani:
            results['ani'] = ani
        else:
            notes.append('ANI skipped: neither fastANI nor skani is installed.')
        return results, notes

    @staticmethod
    def _fetch_genome(accession, work):
        """Download one assembly with NCBI datasets; returns the fasta path or ''."""
        existing = glob.glob(os.path.join(work, accession, '**', '*.fna'), recursive=True)
        if existing:
            return existing[0]
        archive = os.path.join(work, f'{accession}.zip')
        done = _run(['datasets', 'download', 'genome', 'accession', accession,
                     '--include', 'genome', '--filename', archive])
        if done.returncode != 0 or not os.path.exists(archive):
            return ''
        try:
            with zipfile.ZipFile(archive) as z:
                z.extractall(os.path.join(work, accession))
        except zipfile.BadZipFile:
            return ''
        finally:
            if os.path.exists(archive):
                os.remove(archive)
        fastas = glob.glob(os.path.join(work, accession, '**', '*.fna'), recursive=True)
        return fastas[0] if fastas else ''

    @staticmethod
    def _ani(query, reference, work, cpu):
        """Average nucleotide identity between two assemblies, via fastANI or skani."""
        if shutil.which('fastANI'):
            out = os.path.join(work, 'fastani.txt')
            _run(['fastANI', '-q', query, '-r', reference, '-o', out, '-t', str(cpu)])
            if os.path.exists(out):
                with open(out) as f:
                    fields = f.readline().split('\t')
                if len(fields) >= 3:
                    return f'{float(fields[2]):.2f}%'
        if shutil.which('skani'):
            done = subprocess.run(['skani', 'dist', '-t', str(cpu), query, reference],
                                  capture_output=True, text=True)
            lines = done.stdout.splitlines()
            if len(lines) > 1:
                fields = lines[1].split('\t')
                if len(fields) >= 3:
                    return f'{float(fields[2]):.2f}%'
        return ''

    # ------------------------------------------------------------------ MLST
    @staticmethod
    def mlst(sample, assembly):
        """Classic 7-gene typing; returns ({'scheme':, 'st':, 'alleles':}, notes)."""
        if _missing('mlst'):
            return {}, ['MLST skipped: not installed (conda install mlst).']
        done = subprocess.run(['mlst', '--nopath', '--quiet', assembly], capture_output=True, text=True)
        fields = done.stdout.strip().split('\t')
        if done.returncode != 0 or len(fields) < 3:
            return {}, ['MLST failed: no scheme matched.']
        return {'scheme': fields[1], 'st': fields[2], 'alleles': ' '.join(fields[3:])}, []

    # ------------------------------------------------------------------ MOB-suite
    @staticmethod
    def mob_suite(sample, assembly, work, cpu):
        """Plasmid reconstruction and typing; returns (one row per plasmid, notes)."""
        if _missing('mob_recon'):
            return [], ['Plasmid typing skipped: MOB-suite is not installed '
                        '(conda install mob_suite, then hashbrown db mob).']
        out = os.path.join(work, 'mob_recon')
        results = os.path.join(out, 'mobtyper_results.txt')
        contig_report = os.path.join(out, 'contig_report.txt')
        if not os.path.exists(contig_report):
            _run(['mob_recon', '-i', assembly, '-o', out, '-n', str(cpu), '--force'])
        if not os.path.exists(contig_report):
            return [], ['Plasmid typing failed: see MOB-suite logs '
                        '(its database comes with hashbrown db mob).']
        if not os.path.exists(results):
            return [], ['MOB-suite: no plasmids detected.']
        rows = []
        with open(results) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                rows.append({'plasmid': row.get('sample_id', '-').split(':')[-1],
                             'size': row.get('size', '-'),
                             'rep': row.get('rep_type(s)', '-') or '-',
                             'relaxase': row.get('relaxase_type(s)', '-') or '-',
                             'mobility': row.get('predicted_mobility', '-')})
        return rows, []

    # ------------------------------------------------------------------ AMR + virulence
    @staticmethod
    def amr_virulence(sample, assembly, work, cpu):
        """AMRFinderPlus plus abricate (VFDB, PlasmidFinder), merged into one gene table."""
        rows, notes = [], []

        def field(row, *names):
            for name in names:
                if row.get(name) not in (None, ''):
                    return row[name]
            return '-'

        if _missing('amrfinder'):
            notes.append('AMR genes skipped: AMRFinderPlus is not installed '
                         '(conda install ncbi-amrfinderplus, then hashbrown db amrfinder).')
        else:
            out = os.path.join(work, 'amrfinder.tsv')
            if not os.path.exists(out):
                _run(['amrfinder', '-n', assembly, '--plus', '--threads', str(cpu), '-o', out])
            if not os.path.exists(out):
                notes.append('AMRFinderPlus failed: is its database installed? '
                             '(hashbrown db amrfinder)')
            else:
                with open(out) as f:
                    for row in csv.DictReader(f, delimiter='\t'):
                        rows.append({
                            'gene': field(row, 'Gene symbol', 'Element symbol'),
                            'source': 'AMRFinderPlus',
                            'identity': field(row, '% Identity to reference sequence',
                                              '% Identity to reference'),
                            'coverage': field(row, '% Coverage of reference sequence',
                                              '% Coverage of reference'),
                            'contig': field(row, 'Contig id'),
                            'info': ' / '.join(x for x in (field(row, 'Class'), field(row, 'Subclass'))
                                               if x != '-') or field(row, 'Sequence name'),
                        })

        if _missing('abricate'):
            notes.append('Virulence/replicon screens skipped: abricate is not installed '
                         '(conda install abricate).')
        else:
            for db, label in (('vfdb', 'VFDB'), ('plasmidfinder', 'PlasmidFinder')):
                out = os.path.join(work, f'abricate_{db}.tsv')
                if not os.path.exists(out):
                    with open(out + '.part', 'w') as f:
                        done = subprocess.run(['abricate', '--db', db, '--nopath', assembly],
                                              stdout=f, stderr=subprocess.DEVNULL)
                    if done.returncode == 0:
                        os.rename(out + '.part', out)
                if not os.path.exists(out):
                    notes.append(f'abricate --db {db} failed.')
                    continue
                with open(out) as f:
                    text = f.read().lstrip('#')
                for row in csv.DictReader(text.splitlines(), delimiter='\t'):
                    rows.append({'gene': field(row, 'GENE'), 'source': label,
                                 'identity': field(row, '%IDENTITY'),
                                 'coverage': field(row, '%COVERAGE'),
                                 'contig': field(row, 'SEQUENCE'),
                                 'info': field(row, 'PRODUCT', 'RESISTANCE')})
        return rows, notes

    # ------------------------------------------------------------------ annotation
    @staticmethod
    def annotate(sample, assembly, work, cpu, mode='auto', organism=''):
        """Genome annotation: Bakta by default, PGAP on request; returns (result dict, notes)."""
        if mode == 'none':
            return {}, []
        if mode == 'pgap':
            return ExtraQcMethods._pgap(sample, assembly, work, cpu, organism)
        if _missing('bakta'):
            return {}, ['Annotation skipped: Bakta is not installed '
                        '(conda install bakta, then hashbrown db bakta).']
        db = os.environ.get('BAKTA_DB', '')
        if not db:
            versions = glob.glob(os.path.join(db_dir(), 'bakta', '**', 'version.json'),
                                 recursive=True)
            db = os.path.dirname(versions[0]) if versions else ''
        if not db:
            return {}, ['Annotation skipped: no Bakta database (run hashbrown db bakta).']
        out = os.path.join(work, 'bakta')
        summary_file = os.path.join(out, f'{sample}.txt')
        if not os.path.exists(summary_file):
            _run(['bakta', '--db', db, '--output', out, '--prefix', sample,
                  '--threads', str(cpu), '--force', '--skip-plot', assembly])
        if not os.path.exists(summary_file):
            return {}, ['Bakta failed: see its logs.']
        counts = {}
        with open(summary_file) as f:
            for line in f:
                match = re.match(r'^(CDSs|tRNAs|rRNAs|ncRNAs|CRISPR arrays):\s*(\d+)', line.strip())
                if match:
                    counts[match.group(1)] = match.group(2)
        summary = ', '.join(f'{count} {name}' for name, count in counts.items()) or 'done'
        gff = os.path.join(out, f'{sample}.gff3')
        return {'tool': 'Bakta', 'summary': summary,
                'gff': gff if os.path.exists(gff) else ''}, []

    @staticmethod
    def _pgap(sample, assembly, work, cpu, organism):
        """NCBI PGAP via its pgap.py wrapper (Docker); submission-grade but slow. Experimental."""
        if _missing('pgap.py'):
            return {}, ['PGAP skipped: pgap.py is not on the PATH '
                        '(github.com/ncbi/pgap; it needs Docker and ~30 GB).']
        if not organism or organism == '-':
            return {}, ['PGAP skipped: it needs the species name — '
                        'enable the species ID module (hashbrown db mash) first.']
        out = os.path.join(work, 'pgap')
        gff = os.path.join(out, 'annot.gff')
        if not os.path.exists(gff):
            _run(['pgap.py', '--no-self-update', '--cpus', str(cpu), '-o', out,
                  '-g', assembly, '-s', organism])
        if not os.path.exists(gff):
            return {}, ['PGAP failed: see its logs in the report qc folder.']
        return {'tool': 'PGAP', 'summary': 'finished (annot.gff)', 'gff': gff}, []
