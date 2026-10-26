"""Assembly QC: polishing dot plots (nucmer), mapping QC (Qualimap), metrics (QUAST), coverage."""

import csv
import glob
import gzip
import os
import shutil
import subprocess
from concurrent import futures

from hashbrown.methods import Methods


class AssemblyQcMethods:
    @staticmethod
    def fasta_length(input_fasta):
        if not input_fasta or not os.path.exists(input_fasta):
            return 0

        opener = gzip.open if input_fasta.endswith('.gz') else open
        with opener(input_fasta, 'rt') as f:
            return sum(len(line.strip()) for line in f if not line.startswith('>'))

    @staticmethod
    def run_nucmer(sample, ref, query, output_folder, suffix, cpu):
        """Align two assemblies of the same sample and draw a dot plot of the differences."""
        if not (ref and query and os.path.exists(ref) and os.path.exists(query)):
            return

        output_name = os.path.join(output_folder, f'{sample}_{suffix}')

        cmd_nucmer = ['nucmer',
                      '-t', str(cpu),
                      '-p', output_name,  # appends '.delta'
                      '--mincluster=100',
                      ref, query]
        cmd_plot = ['mummerplot',
                    '-x', f'[0,{AssemblyQcMethods.fasta_length(ref)}]',
                    '-y', f'[0,{AssemblyQcMethods.fasta_length(query)}]',
                    '-postscript',
                    '-p', output_name,
                    output_name + '.delta']

        subprocess.run(cmd_nucmer, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        Methods.fix_mummerplot()  # the gnuplot path is undefined in conda's mummerplot
        subprocess.run(cmd_plot, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_nucmer_parallel(sample_dict, output_folder, stage, cpu, parallel):
        """Dot plots before/after polishing: stage 'medaka' (raw vs medaka) or 'polypolish' (medaka vs final)."""
        Methods.make_folder(output_folder)

        def unoriented(path):
            # On a resumed run the final assembly is already dnaapler-reoriented; compare the
            # pre-reorientation sidecar instead, or the rotation would look like a structural change
            sidecar = path.replace('.fasta', '_unoriented.fasta') if path else path
            return sidecar if path and os.path.exists(sidecar) else path

        def pair(sample_obj):
            if stage == 'medaka':
                return unoriented(sample_obj.assembly.raw), unoriented(sample_obj.assembly.medaka)
            return (unoriented(sample_obj.assembly.medaka) or unoriented(sample_obj.assembly.raw),
                    unoriented(sample_obj.assembly.polypolish))

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, *pair(sample_obj), output_folder, stage, max(1, int(cpu / parallel)))
                    for sample, sample_obj in sample_dict.items())
            for _ in executor.map(lambda x: AssemblyQcMethods.run_nucmer(*x), args):
                pass

        # Keep the plots, drop the alignment working files
        for ext in ('.delta', '.fplot', '.gp', '.rplot'):
            for working_file in glob.glob(os.path.join(output_folder, '*' + ext)):
                os.remove(working_file)

    @staticmethod
    def map_minimap2(reads, read_type, ref, output_folder, prefix, cpu):
        """Map one stage's reads onto the assembly; returns the sorted, deduplicated bam path."""
        reads = [f for f in reads if f and os.path.exists(f)]
        if not (reads and ref and os.path.exists(ref)):
            return ''
        r1 = reads[0]
        r2 = reads[1] if len(reads) == 2 else ''

        output_bam = os.path.join(output_folder, prefix + '.bam')

        preset = 'map-ont' if read_type == 'nanopore' else 'sr'
        minimap2_cmd = ['minimap2', '-a', '-x', preset, '-t', str(cpu), ref, r1] + ([r2] if r2 else [])
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-F', '4', '-h', '-']
        samtools_fixmate_cmd = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '-']
        samtools_markdup_cmd = ['samtools', 'markdup', '-r', '-@', str(cpu), '-', output_bam]
        samtools_index_cmd = ['samtools', 'index', output_bam]

        # minimap2 | view (drop unmapped) | [fixmate (pairs only)] | sort | markdup
        steps = [minimap2_cmd, samtools_view_cmd]
        if r2:
            steps.append(samtools_fixmate_cmd)
        steps += [samtools_sort_cmd, samtools_markdup_cmd]

        processes = list()
        previous_stdout = None
        for step in steps:
            p = subprocess.Popen(step, stdin=previous_stdout, stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL)
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = p.stdout
            processes.append(p)
        processes[-1].communicate()

        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        if not os.path.exists(output_bam):
            output_bam = ''

        return output_bam

    @staticmethod
    def mapping_qc(sample, sample_obj, read_type, output_folder, cpu):
        """flagstat + mosdepth for EVERY read-processing stage present (raw, trimmed, filtered).

        Writes <sample>_<stage>.flagstat.txt and <sample>_<stage>.mosdepth.summary.txt per stage,
        and keeps the most processed stage under the legacy names <sample>.flagstat.txt /
        <sample>.mosdepth.summary.txt (what the per-contig coverage uses). Returns
        (sample, {stage: mean depth}).
        """
        if read_type == 'nanopore':
            stage_reads = [('raw', [sample_obj.nanopore.raw]),
                           ('trimmed', [sample_obj.nanopore.trimmed]),
                           ('filtered', [sample_obj.nanopore.filtered])]
        else:  # illumina
            stage_reads = [('raw', list(sample_obj.illumina.raw)),
                           ('trimmed', list(sample_obj.illumina.trimmed))]
        ref = sample_obj.best_assembly()

        depths = {}
        best_prefix = ''
        for stage, reads in stage_reads:
            prefix = os.path.join(output_folder, f'{sample}_{stage}')
            bam = AssemblyQcMethods.map_minimap2(reads, read_type, ref, output_folder,
                                                 f'{sample}_{stage}', cpu)
            if not bam:
                continue

            with open(prefix + '.flagstat.txt', 'w') as f:
                subprocess.run(['samtools', 'flagstat', bam], stdout=f, stderr=subprocess.DEVNULL)
            subprocess.run(['mosdepth', '-x', '-n', '-t', str(cpu), prefix, bam],
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # The overall mean is on the "total" row of the summary
            summary = prefix + '.mosdepth.summary.txt'
            if os.path.exists(summary):
                with open(summary) as f:
                    for line in f:
                        fields = line.split('\t')
                        if fields[0] == 'total':
                            depths[stage] = round(float(fields[3]))
                            best_prefix = prefix

            # The bam and mosdepth working files were only needed here
            for ext in ('.bam', '.bam.bai', '.bai', '.mosdepth.global.dist.txt'):
                leftover = prefix + ext
                if os.path.exists(leftover):
                    os.remove(leftover)

        # The most processed stage also under the legacy names
        if best_prefix:
            for ext in ('.flagstat.txt', '.mosdepth.summary.txt'):
                if os.path.exists(best_prefix + ext):
                    shutil.copyfile(best_prefix + ext, os.path.join(output_folder, sample + ext))

        return sample, depths

    @staticmethod
    def mapping_qc_parallel(sample_dict, output_folder, read_type, cpu, parallel):
        """Per sample: flagstat + mosdepth reports, plus one coverage table for the run."""
        Methods.make_folder(output_folder)
        cov_file = os.path.join(output_folder,
                                ('long' if read_type == 'nanopore' else 'short') + '_read_cov.tsv')

        rows = list()
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, read_type, output_folder, max(1, int(cpu / parallel)))
                    for sample, sample_obj in sample_dict.items())
            for sample, depths in executor.map(lambda x: AssemblyQcMethods.mapping_qc(*x), args):
                if not depths:
                    rows.append({'Sample': sample, 'Stage': '-', 'AverageCoverage': 'NA'})
                for stage, avg_cov in depths.items():
                    rows.append({'Sample': sample, 'Stage': stage, 'AverageCoverage': avg_cov})

        with open(cov_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Sample', 'Stage', 'AverageCoverage'], delimiter='\t')
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: (r['Sample'], r['Stage'])))

    @staticmethod
    def run_quast(sample_dict, output_folder, cpu):
        assembly_list = [sample_obj.best_assembly() for sample_obj in sample_dict.values()
                         if sample_obj.best_assembly() and os.path.exists(sample_obj.best_assembly())]
        if not assembly_list:
            print('\tNo assembly to run QUAST on.')
            return

        cmd = ['quast.py',
               '--output-dir', output_folder,
               '--threads', str(cpu),
               '--min-contig', str(1000),
               '--no-icarus']
        cmd += assembly_list
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
