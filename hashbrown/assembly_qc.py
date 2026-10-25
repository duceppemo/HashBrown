"""Assembly QC: polishing dot plots (nucmer), mapping QC (Qualimap), metrics (QUAST), coverage."""

import csv
import glob
import gzip
import os
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

        def pair(sample_obj):
            if stage == 'medaka':
                return sample_obj.assembly.raw, sample_obj.assembly.medaka
            return (sample_obj.assembly.medaka or sample_obj.assembly.raw), sample_obj.assembly.polypolish

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
    def map_minimap2(sample, sample_obj, read_type, output_folder, cpu):
        """Map the best reads of that type onto the best assembly; returns the sorted, deduplicated bam."""
        # Figure out which reads to use
        r1, r2 = '', ''
        if read_type == 'nanopore':
            r1 = sample_obj.best_nanopore()
        else:  # illumina
            reads = sample_obj.best_illumina()
            if not reads:
                return sample, ''
            r1 = reads[0]
            r2 = reads[1] if len(reads) == 2 else ''

        ref = sample_obj.best_assembly()
        if not (r1 and ref and os.path.exists(ref)):
            return sample, ''

        output_bam = os.path.join(output_folder, sample + '.bam')

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

        return sample, output_bam

    @staticmethod
    def mapping_qc(sample, sample_obj, read_type, output_folder, cpu):
        """Map the reads once, then flagstat (mapping QC) and mosdepth (depth); returns the mean depth."""
        sample, bam = AssemblyQcMethods.map_minimap2(sample, sample_obj, read_type, output_folder, cpu)
        if not bam:
            return sample, 'NA'

        # Mapping QC
        with open(os.path.join(output_folder, sample + '.flagstat.txt'), 'w') as f:
            subprocess.run(['samtools', 'flagstat', bam], stdout=f, stderr=subprocess.DEVNULL)

        # Depth per contig and overall (writes <sample>.mosdepth.summary.txt)
        prefix = os.path.join(output_folder, sample)
        cmd_mosdepth = ['mosdepth', '-x', '-n', '-t', str(cpu), prefix, bam]
        subprocess.run(cmd_mosdepth, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # The overall mean is on the "total" row of the summary
        avg_cov = 0
        summary = prefix + '.mosdepth.summary.txt'
        if os.path.exists(summary):
            with open(summary) as f:
                for line in f:
                    fields = line.split('\t')
                    if fields[0] == 'total':
                        avg_cov = round(float(fields[3]))

        # The bam and mosdepth working files were only needed here
        for ext in ('.bam', '.bam.bai', '.bai', '.mosdepth.global.dist.txt'):
            leftover = os.path.join(output_folder, sample + ext)
            if os.path.exists(leftover):
                os.remove(leftover)

        return sample, avg_cov

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
            for sample, avg_cov in executor.map(lambda x: AssemblyQcMethods.mapping_qc(*x), args):
                rows.append({'Sample': sample, 'AverageCoverage': avg_cov})

        with open(cov_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Sample', 'AverageCoverage'], delimiter='\t')
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: r['Sample']))

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
