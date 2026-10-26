"""Command line interface and pipeline orchestration."""

import os
import sys
from argparse import ArgumentParser
from multiprocessing import cpu_count

from psutil import virtual_memory

from hashbrown import __version__
from hashbrown.assembly_qc import AssemblyQcMethods
from hashbrown.illumina import IlluminaMethods
from hashbrown.methods import Methods
from hashbrown.nanopore import NanoporeMethods
from hashbrown.report import ReportMethods

LONG_ASSEMBLERS = ('flye', 'shasta')
SHORT_ASSEMBLERS = ('skesa', 'spades')
SHORT_POLISHERS = {'polypolish+pypolca': ('polypolish', 'pypolca'),
                   'polypolish': ('polypolish',),
                   'pypolca': ('pypolca',),
                   'none': ()}


def parse_arguments(argv=None):
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1e9)  # in GB

    parser = ArgumentParser(prog='hashbrown',
                            description='HashBrown: long-read, short-read and hybrid bacterial genome '
                                        'assembly pipeline.')
    parser.add_argument('-l', '--long-reads', metavar='/path/to/nanopore_folder',
                        type=str,
                        help='Folder that contains the Nanopore fastq files, one sample per file. '
                             'Required unless assembling short reads only.')
    parser.add_argument('-s', '--short-reads', metavar='/path/to/illumina_folder',
                        type=str,
                        help='Folder that contains the Illumina fastq files ("_R1"/"_R2" mark the pairs; '
                             'a file without them is single-end). With long reads, the short reads polish '
                             'the assemblies and sample names must match (everything before the first "_"). '
                             'Alone, the samples are assembled from the short reads.')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/',
                        help='Folder to hold the result files. Required.')
    parser.add_argument('-a', '--assembler',
                        choices=LONG_ASSEMBLERS + SHORT_ASSEMBLERS,
                        type=str,
                        help='Assembler. With long reads: "flye" (default) or "shasta". '
                             'Short reads only: "skesa" (default) or "spades".')
    parser.add_argument('--min-size', metavar='3000',
                        type=int,
                        help='Minimum read length for Shasta (default 3,000) or minimum read overlap '
                             'for Flye (default "auto").')
    parser.add_argument('--size', metavar='5000000',
                        type=int,
                        help='Override the automatically detected genome size for Flye. Also used by the '
                             'long read filtering step to keep the best 100x, when coverage allows.')
    parser.add_argument('--trim-long',
                        action='store_true',
                        help='Trim long reads with Porechop before assembling. Default is off.')
    parser.add_argument('--filter-long',
                        action='store_true',
                        help='Filter long reads with Filtlong before assembling (drops the worst 5%%). '
                             'Default is off.')
    parser.add_argument('--trim-short',
                        action='store_true',
                        help='Trim short reads with fastp before polishing or assembling. Default is off.')
    parser.add_argument('--long-polisher',
                        choices=['medaka', 'none'], default='medaka',
                        help='Long read polisher for the long read assembly. Default is "medaka". '
                             'Use "none" for non-Nanopore long reads: Medaka is Nanopore-trained and '
                             'actively damages assemblies made from other platforms.')
    parser.add_argument('--model',
                        type=str, default='r1041_e82_400bps_bacterial_methylation',
                        help='Medaka model. The default is ONT\'s bacterial model (R10.4.1 flowcells, '
                             'any recent basecaller), the best performer in our benchmarks. For other '
                             'flowcells pick the model matching your basecaller (--list-models).')
    parser.add_argument('--list-models',
                        action='store_true',
                        help='Show the available Medaka models and exit.')
    parser.add_argument('--short-polisher',
                        choices=list(SHORT_POLISHERS), default='polypolish+pypolca',
                        help='Short read polisher(s) for the long read assembly, applied in order. '
                             'Default is "%(default)s" (pyPOLCA runs with --careful).')
    parser.add_argument('--ion-torrent',
                        action='store_true',
                        help='The short reads are Ion Torrent data (SPAdes assemblies only).')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        type=int, default=max_cpu,
                        help=f'Number of threads. Default is all ({max_cpu}).')
    parser.add_argument('-p', '--parallel', metavar='2',
                        type=int, default=2,
                        help='Number of samples to process in parallel. Keep low on a low-memory '
                             'computer. Default is 2.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        type=int, default=max_mem,
                        help=f'Memory in GB. Default is 85%% of the total memory ({max_mem}).')
    parser.add_argument('-v', '--version', action='version',
                        version=f'HashBrown v{__version__}')

    args = parser.parse_args(argv)

    if args.list_models:
        return args

    # Which pipeline, and is the flag combination coherent?
    if not args.output:
        parser.error('the output folder (-o) is required.')
    if not args.long_reads and not args.short_reads:
        parser.error('provide long reads (-l), short reads (-s), or both.')
    if not args.assembler:
        args.assembler = 'flye' if args.long_reads else 'skesa'
    if args.long_reads and args.assembler in SHORT_ASSEMBLERS:
        parser.error(f'{args.assembler} assembles short reads only: with long reads, '
                     'use "flye" or "shasta" (the short reads are used for polishing).')
    if not args.long_reads and args.assembler in LONG_ASSEMBLERS:
        parser.error(f'{args.assembler} needs long reads (-l).')
    if not args.long_reads and (args.trim_long or args.filter_long or args.min_size or args.size):
        parser.error('--trim-long, --filter-long, --min-size and --size need long reads (-l).')
    if args.ion_torrent and args.assembler != 'spades':
        parser.error('--ion-torrent is only used by SPAdes (-a spades).')

    return args


class HashBrown:
    def __init__(self, args):
        # I/O
        self.long_reads = os.path.abspath(args.long_reads) if args.long_reads else ''
        self.short_reads = os.path.abspath(args.short_reads) if args.short_reads else ''
        self.output_folder = os.path.abspath(args.output)

        # Options
        self.assembler = args.assembler
        self.ref_size = args.size
        self.min_size = args.min_size
        self.trim_long = args.trim_long
        self.filter_long = args.filter_long
        self.trim_short = args.trim_short
        self.long_polisher = args.long_polisher
        self.model = args.model
        self.short_polishers = SHORT_POLISHERS[args.short_polisher]
        self.ion_torrent = args.ion_torrent

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Data
        self.sample_dict = dict()

    def checks(self):
        print('Checking a few things...')

        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        if self.long_reads:
            Methods.check_input_folder(self.long_reads)
        if self.short_reads:
            Methods.check_input_folder(self.short_reads)

        steps = [self.assembler]
        if self.long_reads:
            steps.append('long')
            if self.long_polisher == 'medaka':
                steps.append('medaka')
            if self.trim_long:
                steps.append('trim_long')
            if self.filter_long:
                steps.append('filter_long')
            if self.short_reads and self.short_polishers:
                steps.append('polish')
        if self.short_reads and self.trim_short:
            steps.append('trim_short')
        Methods.check_dependencies(steps)

        Methods.make_folder(self.output_folder)
        Methods.log_versions(os.path.join(self.output_folder, 'log.txt'))
        print('\tAll good!')

    def run(self):
        self.checks()
        if self.long_reads:
            self.run_long()
        else:
            self.run_short()
        print('Done!')

    def run_long(self):
        """Long read assembly, optionally polished with Medaka and with the short reads."""
        out = self.output_folder

        # Step completion flags, so an interrupted run resumes where it stopped
        done_trimming_long = os.path.join(out, 'done_trimming_long')
        done_filtering = os.path.join(out, 'done_filtering_long')
        done_assembling = os.path.join(out, 'done_assembling_long')
        done_polishing_long = os.path.join(out, 'done_polishing_long')
        done_trimming_short = os.path.join(out, 'done_trimming_short')
        done_polishing_short = os.path.join(out, 'done_polishing_short')

        # Output folders
        trimmed_folder = os.path.join(out, '1_trimmed_long')
        filtered_folder = os.path.join(out, '2_filtered_long')
        assembled_folder = os.path.join(out, '3_assembled_long')
        polished_long_folder = os.path.join(out, '4_polished_long')
        illumina_trimmed_folder = os.path.join(out, '5a_trimmed_short')
        illumina_trimmed_report_folder = os.path.join(out, '5b_trimmed_short_report')
        polished_short_folder = os.path.join(out, '6_polished_short')
        assembly_qc_folder = os.path.join(out, '7_assembly_qc')
        gfa_folder = os.path.join(assembly_qc_folder, 'assembly_graphs')

        # Get the input files
        self.sample_dict = Methods.get_nanopore_files(self.long_reads)

        # Drop the unclassified reads if present: we don't know which barcode they come from
        self.sample_dict = {sample: sample_obj for sample, sample_obj in self.sample_dict.items()
                            if 'unclassified' not in sample}

        polishing = bool(self.short_reads and self.short_polishers)
        if self.short_reads:
            self.sample_dict = Methods.get_illumina_files(self.short_reads, self.sample_dict)

        # Trim
        if self.trim_long:
            print('Trimming long reads with Porechop...')
            NanoporeMethods.run_porechop_parallel(self.sample_dict, trimmed_folder,
                                                  self.cpu, self.parallel, done_trimming_long)

        # Filter
        if self.filter_long:
            print('Filtering long reads with Filtlong...')
            NanoporeMethods.run_filtlong_parallel(self.sample_dict, filtered_folder,
                                                  self.ref_size, self.parallel, done_filtering)

        # Assemble
        print(f'Assembling long reads with {self.assembler.capitalize()}...')
        if self.assembler == 'flye':
            NanoporeMethods.assemble_flye_parallel(self.sample_dict, assembled_folder, gfa_folder,
                                                   self.ref_size, self.min_size, self.cpu, self.parallel,
                                                   done_assembling, out)
        else:  # shasta
            NanoporeMethods.assemble_shasta_parallel(self.sample_dict, assembled_folder, gfa_folder,
                                                     self.min_size, self.cpu, self.parallel,
                                                     done_assembling, out)

        # Long read polishing
        if self.long_polisher == 'medaka':
            print('Long read polishing with Medaka...')
            NanoporeMethods.polish_medaka_loop(self.sample_dict, polished_long_folder, self.model,
                                               self.cpu, done_polishing_long)

        # Short read polishing
        if polishing:
            if self.trim_short:
                print('Trimming short reads with fastp...')
                IlluminaMethods.trim_fastp_parallel(self.sample_dict, illumina_trimmed_folder,
                                                    illumina_trimmed_report_folder,
                                                    self.cpu, self.parallel, done_trimming_short)

            print('Short read polishing with {}...'.format(' and '.join(self.short_polishers)))
            IlluminaMethods.polish_parallel(self.sample_dict, polished_short_folder,
                                            self.short_polishers, self.cpu, self.parallel,
                                            done_polishing_short)

        # QC
        print('Performing assembly QC...')
        quast_folder = os.path.join(assembly_qc_folder, 'quast')
        mapping_qc_long_folder = os.path.join(assembly_qc_folder, 'mapping_qc_long_reads')
        mapping_qc_short_folder = os.path.join(assembly_qc_folder, 'mapping_qc_short_reads')
        mummer_long_folder = os.path.join(assembly_qc_folder, 'mummer_long')
        mummer_short_folder = os.path.join(assembly_qc_folder, 'mummer_short')

        if self.long_polisher == 'medaka':
            print('\tDot plots - raw vs long read polished (nucmer)')
            AssemblyQcMethods.run_nucmer_parallel(self.sample_dict, mummer_long_folder, 'medaka',
                                                  self.cpu, self.parallel)

        print('\tQUAST')
        AssemblyQcMethods.run_quast(self.sample_dict, quast_folder, self.cpu)

        print('\tMapping QC and coverage - long reads (flagstat, mosdepth)')
        AssemblyQcMethods.mapping_qc_parallel(self.sample_dict, mapping_qc_long_folder, 'nanopore',
                                              self.cpu, self.parallel)

        if polishing:
            print('\tDot plots - long read polished vs short read polished (nucmer)')
            AssemblyQcMethods.run_nucmer_parallel(self.sample_dict, mummer_short_folder, 'polypolish',
                                                  self.cpu, self.parallel)

            print('\tMapping QC and coverage - short reads (flagstat, mosdepth)')
            AssemblyQcMethods.mapping_qc_parallel(self.sample_dict, mapping_qc_short_folder, 'illumina',
                                                  self.cpu, self.parallel)

        # Report
        print('Building the report...')
        ReportMethods.build_reports(self.sample_dict, out, os.path.join(out, '8_report'), 'long', self.cpu)

    def run_short(self):
        """Short read assembly with SKESA or SPAdes (single-end or paired-end)."""
        out = self.output_folder

        done_trimming_short = os.path.join(out, 'done_trimming_short')
        done_assembling = os.path.join(out, 'done_assembling_short')

        trimmed_folder = os.path.join(out, '1_trimmed_short')
        trimmed_report_folder = os.path.join(out, '1_trimmed_short', 'report')
        assembled_folder = os.path.join(out, '2_assembled_short')
        assembly_qc_folder = os.path.join(out, '3_assembly_qc')
        gfa_folder = os.path.join(assembly_qc_folder, 'assembly_graphs')

        # Get the input files
        self.sample_dict = Methods.get_illumina_files(self.short_reads, self.sample_dict,
                                                      require_long_reads=False)

        # Trim
        if self.trim_short:
            print('Trimming short reads with fastp...')
            IlluminaMethods.trim_fastp_parallel(self.sample_dict, trimmed_folder, trimmed_report_folder,
                                                self.cpu, self.parallel, done_trimming_short)

        # Assemble
        assembler_name = 'SKESA' if self.assembler == 'skesa' else 'SPAdes'
        print(f'Assembling short reads with {assembler_name}...')
        IlluminaMethods.assemble_short_parallel(self.sample_dict, self.assembler, assembled_folder,
                                                gfa_folder, self.ion_torrent, self.cpu, self.mem,
                                                self.parallel, done_assembling)

        # QC
        print('Performing assembly QC...')
        quast_folder = os.path.join(assembly_qc_folder, 'quast')
        mapping_qc_short_folder = os.path.join(assembly_qc_folder, 'mapping_qc_short_reads')

        print('\tQUAST')
        AssemblyQcMethods.run_quast(self.sample_dict, quast_folder, self.cpu)

        print('\tMapping QC and coverage - short reads (flagstat, mosdepth)')
        AssemblyQcMethods.mapping_qc_parallel(self.sample_dict, mapping_qc_short_folder, 'illumina',
                                              self.cpu, self.parallel)

        # Report
        print('Building the report...')
        ReportMethods.build_reports(self.sample_dict, out, os.path.join(out, '4_report'), 'short', self.cpu)


def main(argv=None):
    args = parse_arguments(argv)

    if args.list_models:
        models = NanoporeMethods.list_medaka_models()
        if models:
            print('Available Medaka models:\n{}'.format(', '.join(models)))
        else:
            print('Could not list the Medaka models. Is medaka installed?')
        return

    try:
        HashBrown(args).run()
    except Exception as error:
        sys.exit(f'Error: {error}')


if __name__ == '__main__':
    main()
