import os
import shutil
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from bga_short_methods import Methods


__author__ = 'duceppemo'
__version__ = '0.1'


'''
mamba create -n short_BGA falco=1.2.1 skesa=2.4 spades=3.15.5 bbmap=39.01 fastp=0.23.3 samtools=1.17 minimap2=2.26 \
    psutil=5.9.5 qualimap=2.2.2d quast=5.2.0 bandage=0.8.1
'''


class BGAShort(object):
    def __init__(self, args):
        # I/O
        self.input = args.input
        self.output = args.output

        # Args
        self.assembler = args.assembler

        # Performance
        self.cpu = args.threads
        self.mem = args.memory
        self.parallel = args.parallel

        # Script location
        script_home = __file__
        self.gfa_con = '/'.join(script_home.split('/')[:-1]) + '/bin/gfa_connector'

        # Run pipeline
        self.run()

    def run(self):
        # Checks
        print('Checking a few things...')
        Methods.check_cpus(self.cpu, self.parallel)
        Methods.check_mem(self.mem)
        Methods.check_input(self.input)

        # Folder structure
        read_qc_folder = self.output + '/1_read_qc/'
        trimmed_folder = self.output + '/2_trimmed/'
        assembled_folder = self.output + '/3_assembled/'
        assembly_qc_folder = self.output + '/4_assembly_qc/'
        quast_folder = assembly_qc_folder + 'quast/'
        qualimap_folder = assembly_qc_folder + 'qualimap/'
        gfa_folder = assembly_qc_folder + 'assembly_graphs/'

        # Create folders
        Methods.make_folder(read_qc_folder)
        Methods.make_folder(trimmed_folder)
        Methods.make_folder(assembled_folder)
        Methods.make_folder(assembly_qc_folder)
        Methods.make_folder(quast_folder)
        Methods.make_folder(qualimap_folder)
        Methods.make_folder(gfa_folder)

        # Get fastq files
        sample_dict = Methods.get_fastq_files(self.input)

        # Read QC
        print('Performing read QC with FastQC...')
        Methods.read_qc_falco_parallel(sample_dict, read_qc_folder, self.cpu, self.parallel)

        # Trimming
        print('Trimming reads with FastP...')
        Methods.trim_fastp_parallel(sample_dict, trimmed_folder, trimmed_folder + 'report/', self.cpu, self.parallel)

        # Update dictionary
        sample_dict = Methods.get_fastq_files(trimmed_folder)

        # Assemble
        print('Assembling reads with {}...'.format(self.assembler))
        if self.assembler == 'skesa':
            Methods.assemble_skesa_parallel(sample_dict, assembled_folder, self.cpu, self.mem, self.parallel)
        else:  # elif self.assembler == 'spades':
            Methods.assemble_spades_parallel(sample_dict, assembled_folder, self.cpu, self.mem, self.parallel)

        # Assembly QC
        print('Performing assembly QC...')
        assembly_list = Methods.list_files_in_folder(assembled_folder, '.fasta')

        print('\tMapping trimmed reads to assemblies with minimap2')
        Methods.map_minimap2_parallel(assembly_list, trimmed_folder, qualimap_folder, self.cpu, self.parallel)
        bam_list = Methods.list_files_in_folder(qualimap_folder, '.bam')

        print('\tRunning Qualimap')
        Methods.run_qualimap_parallel(bam_list, qualimap_folder, self.cpu, self.mem, self.parallel)

        print('Running QUAST...')
        Methods.run_quast(assembly_list, quast_folder, self.cpu)

        print('Drawing assembly graphs...')
        if self.assembler == 'skesa':
            Methods.gfa_connector_parallel(assembly_list, sample_dict, gfa_folder, self.cpu, self.parallel,
                                           self.gfa_con)
        else:  # elif self.assembler == 'spades':
            # Move gfa files
            gfa_list = Methods.list_files_in_folder(assembled_folder, '.gfa')
            for gfa in gfa_list:
                shutil.move(gfa, gfa_folder)

        gfa_list = Methods.list_files_in_folder(gfa_folder, '.gfa')
        Methods.assembly_graph_parallel(gfa_list, gfa_folder, self.cpu)

        print('DONE')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Bacterial Genome assembly for short reads.')
    parser.add_argument('-i', '--input', metavar='/path/to/short_reads',
                        required=True, type=str,
                        help='Folder that contains the sort read files in fastq format. Gzipped or not. Single-end or '
                             'paired-end. Samples will be named according to everything before the first '
                             'underscore ("_"). Mandatory')
    parser.add_argument('--ion-torrent',
                        required=False, action='store_true',
                        help='Use this flag if fastq files are from Ion Torrent. Needed for SPAdes. Optional')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-a', '--assembler',
                        required=False, default='skesa',
                        choices=['skesa', 'spades'],
                        type=str,
                        help='Assembly method. Default "skesa". Optional.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Keep low if your computer has low memory. '
                             'Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    BGAShort(arguments)
