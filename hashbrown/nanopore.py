"""Long-read steps: trimming (Porechop), filtering (Filtlong), assembly (Flye or Shasta), polishing (Medaka)."""

import contextlib
import csv
import gzip
import os
import shutil
import subprocess
from concurrent import futures
from glob import glob

from hashbrown.methods import Methods


class NanoporeMethods:
    @staticmethod
    def list_medaka_models():
        """The models medaka_consensus accepts (parsed from its help, so always in sync)."""
        try:
            p = subprocess.run(['medaka_consensus', '-h'], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
        except FileNotFoundError:
            return []
        for line in p.stdout.splitlines():
            if 'Choices:' in line:
                return sorted(line.split(':', 1)[1].split())
        return []

    @staticmethod
    def run_porechop(sample, sample_obj, trimmed_folder, cpu, flag):
        input_fastq = sample_obj.nanopore.raw
        trimmed_fastq = os.path.join(trimmed_folder, sample + '.fastq.gz')

        if not os.path.exists(flag):
            cmd = ['porechop',
                   '-i', input_fastq,
                   '-o', trimmed_fastq,
                   '--threads', str(cpu),
                   '--check_reads', str(1000)]

            print(f'\t{sample}')
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # The file may be missing when resuming a partial run
        if not os.path.exists(trimmed_fastq):
            trimmed_fastq = ''

        return sample, trimmed_fastq

    @staticmethod
    def run_porechop_parallel(sample_dict, output_folder, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        # When the flag file exists the trimming is skipped, but the loop still runs
        # to point the samples at the files from the previous run.
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, max(1, int(cpu / parallel)), flag)
                    for sample, sample_obj in sample_dict.items())
            for sample, trimmed in executor.map(lambda x: NanoporeMethods.run_porechop(*x), args):
                sample_dict[sample].nanopore.trimmed = trimmed
                if not trimmed and not os.path.exists(flag):
                    print(f'\tWARNING: Porechop produced nothing for {sample}; '
                          'continuing with the untrimmed reads.')

        if os.path.exists(flag):
            print('\tSkipping trimming long reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def run_filtlong(sample, sample_obj, filtered_folder, genome_size, flag):
        input_fastq = sample_obj.nanopore.trimmed or sample_obj.nanopore.raw
        filtered_fastq = os.path.join(filtered_folder, sample + '.fastq.gz')

        if not os.path.exists(flag):
            cmd = ['filtlong',
                   '--keep_percent', str(95),  # drop the worst 5% of the reads
                   '--min_length', str(500),  # remove rejected reads from targeted sequencing
                   input_fastq]
            if genome_size:
                cmd += ['--target_bases', str(genome_size * 100)]  # keep the best 100x, if more

            # Filtlong writes to stdout
            print(f'\t{sample}')
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            with gzip.open(filtered_fastq, 'wb') as f:  # streamed: the fastq can be many GB
                shutil.copyfileobj(p.stdout, f, length=1024 * 1024)
            p.wait()

        if os.path.exists(filtered_fastq):
            with gzip.open(filtered_fastq, 'rb') as f:
                if not f.read(1):  # filtlong failed silently: an empty file helps nobody
                    os.remove(filtered_fastq)
                    filtered_fastq = ''
        else:
            filtered_fastq = ''

        return sample, filtered_fastq

    @staticmethod
    def run_filtlong_parallel(sample_dict, output_folder, genome_size, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, genome_size, flag)
                    for sample, sample_obj in sample_dict.items())
            for sample, filtered in executor.map(lambda x: NanoporeMethods.run_filtlong(*x), args):
                sample_dict[sample].nanopore.filtered = filtered
                if not filtered and not os.path.exists(flag):
                    print(f'\tWARNING: Filtlong produced nothing for {sample}; '
                          'continuing with the unfiltered reads.')

        if os.path.exists(flag):
            print('\tSkipping filtering long reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def assemble_flye(sample, sample_obj, output_folder, gfa_folder, genome_size, min_size, cpu, flag):
        # I/O
        output_assembly = os.path.join(output_folder, sample + '.fasta')
        output_subfolder = os.path.join(output_folder, sample)  # one working subfolder per sample
        input_fastq = sample_obj.best_nanopore()

        if not os.path.exists(flag):
            cmd_flye = ['flye',
                        '--nano-hq', input_fastq,
                        '--threads', str(cpu),
                        '--out-dir', output_subfolder,
                        '--iterations', str(3)]
            if min_size:
                cmd_flye += ['--min-overlap', str(min_size)]
            if genome_size:
                cmd_flye += ['--genome-size', str(genome_size)]

            print(f'\t{sample}')
            Methods.make_folder(output_subfolder)
            subprocess.run(cmd_flye, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            if os.path.exists(os.path.join(output_subfolder, 'assembly.fasta')):  # assembly succeeded
                shutil.move(os.path.join(output_subfolder, 'assembly.fasta'), output_assembly)

                # Keep the assembly graph, and a picture of it
                Methods.make_folder(gfa_folder)
                gfa_file = os.path.join(gfa_folder, sample + '.gfa')
                shutil.move(os.path.join(output_subfolder, 'assembly_graph.gfa'), gfa_file)
                Methods.make_assembly_graph_image(gfa_file, os.path.join(gfa_folder, sample + '.png'))
            else:
                print(f'\tAssembly failed for {sample}')

        if not os.path.exists(output_assembly):
            output_assembly = ''

        return sample, output_assembly

    @staticmethod
    def assemble_flye_parallel(sample_dict, assembly_folder, gfa_folder, genome_size, min_size, cpu, parallel,
                               flag, output_folder):
        Methods.make_folder(assembly_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, assembly_folder, gfa_folder, genome_size, min_size,
                     max(1, int(cpu / parallel)), flag)
                    for sample, sample_obj in sample_dict.items())
            for sample, assembly in executor.map(lambda x: NanoporeMethods.assemble_flye(*x), args):
                sample_dict[sample].assembly.raw = assembly

        NanoporeMethods.flye_assembly_stats(assembly_folder, output_folder)

        if os.path.exists(flag):
            print('\tSkipping assembling long reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def flye_assembly_stats(assembly_folder, output_folder):
        """Parse each sample's flye.log into one stats table, then drop the Flye working folders."""
        output_stats_file = os.path.join(output_folder, 'flye_stats.tsv')

        log_list = glob(os.path.join(assembly_folder, '**', 'flye.log'), recursive=True)
        if not log_list:  # resuming: the working folders are already cleaned up, keep the existing table
            return

        rows = list()
        for log_file in sorted(log_list):
            sample = os.path.basename(os.path.dirname(log_file))
            stats = {'Sample': sample, 'TotalReadLength': 0, 'ReadsN50': 0,
                     'AssemblyLength': 0, 'Contigs': 0, 'Coverage': 0}
            with open(log_file) as f:
                for line in f:
                    line = line.rstrip()
                    if 'Total read length:' in line:
                        stats['TotalReadLength'] = line.split()[-1]
                    elif 'N50/N90' in line:
                        stats['ReadsN50'] = line.split(':')[-1].split('/')[0].strip()
                    elif 'Total length:' in line:
                        stats['AssemblyLength'] = line.split('\t')[-1]
                    elif 'Fragments:' in line:
                        stats['Contigs'] = line.split('\t')[-1]
                    elif 'Mean coverage:' in line:
                        stats['Coverage'] = line.split('\t')[-1]
            rows.append(stats)

            # The assembly and its graph were moved out already: drop the Flye working folder
            shutil.rmtree(os.path.dirname(log_file), ignore_errors=True)

        with open(output_stats_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def assemble_shasta(sample, sample_obj, output_folder, gfa_folder, min_size, cpu, flag):
        # I/O
        input_fastq = sample_obj.best_nanopore()
        fasta = os.path.join(output_folder, sample + '_reads.fasta')  # Shasta wants uncompressed fasta
        output_assembly = os.path.join(output_folder, sample + '.fasta')
        output_subfolder = os.path.join(output_folder, sample)

        if not os.path.exists(flag):
            cmd_fastq_to_fasta = ['seqtk', 'seq', '-a', input_fastq]
            cmd_shasta = ['shasta',
                          '--config', 'Nanopore-May2022',
                          '--input', fasta,
                          '--assemblyDirectory', output_subfolder,
                          '--threads', str(cpu),
                          '--Reads.minReadLength', str(min_size if min_size else 3000)]
            cmd_shasta_clean = ['shasta',
                                '--assemblyDirectory', output_subfolder,
                                '--command', 'cleanupBinaryData']

            print(f'\t{sample}')

            # Convert fastq to fasta
            with open(fasta, 'w') as f:
                subprocess.run(cmd_fastq_to_fasta, stdout=f)

            # Run the assembler, keeping its stdout: that is where the assembly stats are
            shasta_stdout = os.path.join(output_folder, sample + '_shasta_stdout.txt')
            with open(shasta_stdout, 'w') as f:
                subprocess.run(cmd_shasta, stdout=f, stderr=subprocess.DEVNULL)

            # Cleanup temporary files
            subprocess.run(cmd_shasta_clean, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            os.remove(fasta)

            # The assembly and gfa files are always created, but empty when the assembly failed
            shasta_fasta = os.path.join(output_subfolder, 'Assembly.fasta')
            if os.path.exists(shasta_fasta) and os.stat(shasta_fasta).st_size != 0:
                shutil.move(shasta_fasta, output_assembly)

                Methods.make_folder(gfa_folder)
                gfa_file = os.path.join(gfa_folder, sample + '.gfa')
                shutil.move(os.path.join(output_subfolder, 'Assembly.gfa'), gfa_file)
                Methods.make_assembly_graph_image(gfa_file, os.path.join(gfa_folder, sample + '.png'))
            else:
                print(f'\tAssembly failed for {sample}')

        if not os.path.exists(output_assembly):
            output_assembly = ''

        return sample, output_assembly

    @staticmethod
    def assemble_shasta_parallel(sample_dict, assembly_folder, gfa_folder, min_size, cpu, parallel, flag,
                                 output_folder):
        Methods.make_folder(assembly_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, assembly_folder, gfa_folder, min_size,
                     max(1, int(cpu / parallel)), flag)
                    for sample, sample_obj in sample_dict.items())
            for sample, assembly in executor.map(lambda x: NanoporeMethods.assemble_shasta(*x), args):
                sample_dict[sample].assembly.raw = assembly

        NanoporeMethods.shasta_assembly_stats(assembly_folder, output_folder)

        if os.path.exists(flag):
            print('\tSkipping assembling long reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def shasta_assembly_stats(assembly_folder, output_folder):
        """Parse each sample's Shasta stdout into one stats table, then drop the working files."""
        output_stats_file = os.path.join(output_folder, 'shasta_stats.tsv')

        log_list = glob(os.path.join(assembly_folder, '*_shasta_stdout.txt'))
        if not log_list:  # resuming: keep the existing table
            return

        rows = list()
        for log_file in sorted(log_list):
            sample = os.path.basename(log_file).replace('_shasta_stdout.txt', '')
            stats = {'Sample': sample, 'ReadNumber': 0, 'TotalBases': 0, 'ReadsN50': 0,
                     'AssemblyLength': 0, 'Contigs': 0}
            with open(log_file) as f:
                for line in f:
                    line = line.rstrip()
                    if 'Total number of reads is' in line:
                        stats['ReadNumber'] = line.split()[-1].replace('.', '')
                    elif 'Total number of raw bases is' in line:
                        stats['TotalBases'] = line.split()[-1].replace('.', '')
                    elif 'N50 for read length is' in line:
                        stats['ReadsN50'] = line.split()[-2]
                    elif 'Total length of assembled sequence is' in line:
                        stats['AssemblyLength'] = line.split()[-1]
                    elif 'The assembly graph has' in line:
                        stats['Contigs'] = line.split()[-3]
            rows.append(stats)
            os.remove(log_file)

        with open(output_stats_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

        # Drop the Shasta working folders
        for entry in os.scandir(assembly_folder):
            if entry.is_dir():
                shutil.rmtree(entry.path, ignore_errors=True)

    @staticmethod
    def polish_medaka(sample, sample_obj, output_folder, model, cpu, flag):
        # I/O
        input_fastq = sample_obj.best_nanopore()
        input_assembly = sample_obj.assembly.raw

        sample_subfolder = os.path.join(output_folder, sample)
        default_medaka_output = os.path.join(sample_subfolder, 'consensus.fasta')
        polished_assembly = os.path.join(output_folder, sample + '.fasta')

        if not os.path.exists(flag):
            if input_assembly and os.path.exists(input_assembly):
                Methods.make_folder(sample_subfolder)

                print(f'\t{sample}')
                cmd_medaka = ['medaka_consensus',
                              '-i', input_fastq,
                              '-d', input_assembly,
                              '-o', sample_subfolder,
                              '-t', str(cpu),
                              '-m', model]
                subprocess.run(cmd_medaka, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                if os.path.exists(default_medaka_output):
                    os.rename(default_medaka_output, polished_assembly)

                    # Reformat the fasta to 80 characters per line
                    Methods.format_fasta(polished_assembly, polished_assembly + '.tmp')
                    os.rename(polished_assembly + '.tmp', polished_assembly)

                # Cleanup the working folder and the index files medaka leaves next to the assembly
                shutil.rmtree(sample_subfolder, ignore_errors=True)
                for leftover in (input_assembly + '.map-ont.mmi', input_assembly + '.fai'):
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(leftover)
            else:
                print(f'\tNo assembly for {sample}')

        if not os.path.exists(polished_assembly):
            polished_assembly = ''

        return sample, polished_assembly

    @staticmethod
    def polish_medaka_loop(sample_dict, output_folder, model, cpu, flag):
        """One sample at a time: medaka crashes when several instances run in parallel."""
        Methods.make_folder(output_folder)

        for sample, sample_obj in sample_dict.items():
            sample, polished = NanoporeMethods.polish_medaka(sample, sample_obj, output_folder,
                                                             model, int(cpu), flag)
            sample_dict[sample].assembly.medaka = polished

        if os.path.exists(flag):
            print('\tSkipping long read polishing. Already done.')
        else:
            Methods.flag_done(flag)
