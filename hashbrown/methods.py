"""General helpers: input discovery, resource checks, small file utilities."""

import os
import pathlib
import re
import shutil
import subprocess
import sys
import textwrap
from multiprocessing import cpu_count
from types import SimpleNamespace

from psutil import virtual_memory


class Sample:
    """All the files known for one sample, filled in as the pipeline progresses.

    Empty string (or empty list for paired reads) means "not available (yet)".
    """

    def __init__(self):
        self.nanopore = SimpleNamespace(raw='', trimmed='', filtered='')
        self.illumina = SimpleNamespace(raw=[], trimmed=[])
        self.assembly = SimpleNamespace(raw='', medaka='', polypolish='', unoriented='')
        self.bam = SimpleNamespace(raw='')

    def best_nanopore(self):
        """The most processed long-read file available."""
        return self.nanopore.filtered or self.nanopore.trimmed or self.nanopore.raw

    def best_illumina(self):
        """The most processed short-read files available: [R1, R2], [single] or []."""
        return self.illumina.trimmed or self.illumina.raw

    def best_assembly(self):
        """The most polished assembly available."""
        return self.assembly.polypolish or self.assembly.medaka or self.assembly.raw


class Methods:
    accepted_extensions = ('.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz')

    # External programs each part of the pipeline needs on the PATH
    dependencies = {
        'always': ['samtools', 'minimap2', 'quast.py', 'mosdepth', 'seqkit'],
        'long': ['seqtk', 'nucmer', 'mummerplot', 'gnuplot'],
        'trim_long': ['porechop'],
        'filter_long': ['filtlong'],
        'flye': ['flye'],
        'shasta': ['shasta'],
        'medaka': ['medaka_consensus'],
        'trim_short': ['fastp'],
        'polish': ['bwa', 'polypolish', 'pypolca'],
        'skesa': ['skesa'],
        'spades': ['spades.py'],
    }
    # Missing optional programs only cost their output (a warning says so)
    optional_dependencies = ['Bandage', 'gfa_connector']

    @staticmethod
    def get_sample_name(filename):
        """Everything before the first underscore of the base name (or the whole base name without one)."""
        return os.path.basename(filename).split('.')[0].split('_')[0]

    @staticmethod
    def check_input_folder(my_input):
        if not os.path.isdir(my_input):
            raise Exception(f'"{my_input}" is not an existing folder. '
                            'Please provide a folder containing fastq files.')

        # At least one sequence file must be there (recursively); other files are simply ignored.
        for _root, _, filenames in os.walk(my_input):
            if any(f.endswith(Methods.accepted_extensions) for f in filenames):
                return
        raise Exception(f'No sequence file found in "{my_input}". '
                        f'Accepted extensions: {", ".join(Methods.accepted_extensions)}')

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if not 1 <= requested_cpu <= total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write(f'Number of threads was set to {requested_cpu}\n')
        if not 1 <= n_proc <= requested_cpu:
            n_proc = max(1, min(n_proc, requested_cpu))
            sys.stderr.write(f'Number of samples to process in parallel was set to {n_proc}\n')

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1e9)  # in GB
        if not requested_mem or not 1 <= requested_mem <= max_mem:
            requested_mem = max_mem
            sys.stderr.write(f'Memory was set to {requested_mem}GB\n')
        return requested_mem

    @staticmethod
    def check_dependencies(steps):
        """Make sure every external program the requested steps use is installed."""
        needed = list(Methods.dependencies['always'])
        for step in steps:
            needed += Methods.dependencies[step]

        missing = [tool for tool in needed if not shutil.which(tool)]
        if missing:
            raise Exception('Some required programs are not installed or not on the PATH: {}\n'
                            'Did you activate the HashBrown conda environment?'.format(', '.join(missing)))

        for tool in Methods.optional_dependencies:
            if not shutil.which(tool):
                sys.stderr.write(f'{tool} not found: the corresponding assembly graph output is skipped.\n')

    @staticmethod
    def log_versions(log_file):
        """Best effort: record the version of each installed program for the methods section of a paper."""
        version_commands = {
            'Python': [sys.executable, '--version'],
            'Flye': ['flye', '--version'],
            'Shasta': ['shasta', '--version'],
            'Medaka': ['medaka', '--version'],
            'Porechop': ['porechop', '--version'],
            'Filtlong': ['filtlong', '--version'],
            'Fastp': ['fastp', '--version'],
            'BWA': ['bwa'],
            'Polypolish': ['polypolish', '--version'],
            'pyPOLCA': ['pypolca', '--version'],
            'SAMtools': ['samtools', 'version'],
            'Minimap2': ['minimap2', '--version'],
            'Seqtk': ['seqtk'],
            'QUAST': ['quast.py', '--version'],
            'mosdepth': ['mosdepth', '--version'],
            'SeqKit': ['seqkit', 'version'],
            'MUMmer (nucmer)': ['nucmer', '--version'],
            'Bandage': ['Bandage', '--version'],
        }

        with open(log_file, 'w') as f:
            for name, cmd in version_commands.items():
                if not shutil.which(cmd[0]) and name != 'Python':
                    continue
                try:
                    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       timeout=60, text=True)
                    # The version is on the first line that holds a digit
                    version = next((line.strip() for line in p.stdout.splitlines()
                                    if any(c.isdigit() for c in line)), 'unknown')
                except Exception:
                    version = 'unknown'
                f.write(f'{name}: {version}\n')

    @staticmethod
    def check_read_file(file_path, filename):
        """Fail early and clearly on an unreadable input: a broken symlink or an empty file would
        otherwise ripple through the pipeline as empty intermediate files and a garbage report."""
        if not os.path.isfile(file_path):
            raise Exception(f'input file "{filename}" is a broken link or unreadable ({file_path}).')
        if os.path.getsize(file_path) == 0:
            raise Exception(f'input file "{filename}" is empty ({file_path}).')

    @staticmethod
    def make_folder(folder):
        # Create parent directories if needed; no error if it already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def make_assembly_graph_image(gfa_file, png_file):
        if shutil.which('Bandage'):
            cmd = ['Bandage', 'image', gfa_file, png_file]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def get_nanopore_files(in_folder):
        """Map sample name -> Sample for every sequence file found (recursively) in the folder."""
        sample_dict = dict()

        for root, _, filenames in os.walk(in_folder):
            for filename in sorted(filenames):
                if not filename.endswith(Methods.accepted_extensions):
                    continue
                file_path = os.path.realpath(os.path.join(root, filename))  # follow symbolic links
                Methods.check_read_file(file_path, filename)
                sample = Methods.get_sample_name(filename)

                if sample in sample_dict:
                    sys.stderr.write(f'More than one long-read file for "{sample}": '
                                     f'keeping {sample_dict[sample].nanopore.raw}, ignoring {file_path}\n')
                    continue

                sample_dict[sample] = Sample()
                sample_dict[sample].nanopore.raw = file_path

        if not sample_dict:
            raise Exception(f'No sequence file found in "{in_folder}".')

        return sample_dict

    @staticmethod
    def get_illumina_files(in_folder, sample_dict, require_long_reads=True):
        """Attach the short-read files to the samples.

        "_R1"/"_R2" in the file name mark the two files of a pair; a file without either is a
        single-end sample (both polishers handle single-end reads). With require_long_reads (the
        hybrid pipeline, where the short reads only polish), a short-read sample without matching
        long reads is reported and skipped. Without it (the short-read-only pipeline), every
        sample found is assembled.
        """
        for root, _, filenames in os.walk(in_folder):
            for filename in sorted(filenames):
                if not filename.endswith(Methods.accepted_extensions):
                    continue
                file_path = os.path.realpath(os.path.join(root, filename))
                Methods.check_read_file(file_path, filename)
                sample = Methods.get_sample_name(filename)

                if sample not in sample_dict:
                    if require_long_reads:
                        sys.stderr.write(f'No matching long reads for short-read file {filename}: '
                                         'skipping it.\n')
                        continue
                    sample_dict[sample] = Sample()

                if re.search(r'_R1(?=[._])', filename):  # not _R10: the mate tag ends at . or _
                    sample_dict[sample].illumina.raw.insert(0, file_path)
                else:  # _R2 of a pair, or a single-end sample
                    sample_dict[sample].illumina.raw.append(file_path)

        # A sample is a pair or one single-end file; with more files than that, keep a clean
        # pair if there is exactly one - anything else is ambiguous and must not be guessed at
        for sample, sample_obj in sample_dict.items():
            reads = sample_obj.illumina.raw
            if len(reads) > 2:
                r1 = [f for f in reads if re.search(r'_R1(?=[._])', os.path.basename(f))]
                r2 = [f for f in reads if re.search(r'_R2(?=[._])', os.path.basename(f))]
                if len(r1) == 1 and len(r2) == 1:
                    for extra in reads:
                        if extra not in r1 + r2:
                            sys.stderr.write(f'Too many short-read files for {sample}: '
                                             f'ignoring {os.path.basename(extra)}.\n')
                    sample_obj.illumina.raw = r1 + r2
                else:
                    names = ', '.join(sorted(os.path.basename(f) for f in reads))
                    raise Exception(f'ambiguous short-read files for "{sample}" ({names}): '
                                    'give one _R1/_R2 pair or one single-end file '
                                    '(merge lanes beforehand).')

        if not sample_dict:
            raise Exception(f'No sequence file found in "{in_folder}".')

        return sample_dict

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'w') as f:
            for item in my_list:
                f.write(f'{item}\n')

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w'):
            pass

    @staticmethod
    def format_fasta(input_fasta, output_fasta, width=80):
        """Rewrite a fasta file with sequence lines wrapped at a fixed width."""
        with open(input_fasta) as in_handle, open(output_fasta, 'w') as out_handle:
            seq = list()
            for line in in_handle:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seq:
                        out_handle.write('\n'.join(textwrap.wrap(''.join(seq), width,
                                                                 break_long_words=True)) + '\n')
                        seq = list()
                    out_handle.write(line + '\n')
                else:
                    seq.append(line)
            if seq:
                out_handle.write('\n'.join(textwrap.wrap(''.join(seq), width, break_long_words=True)) + '\n')

    @staticmethod
    def fix_mummerplot():
        """Fix the gnuplot path in mummerplot, left undefined by the conda install."""
        mummerplot_path = shutil.which('mummerplot')
        gnuplot_path = shutil.which('gnuplot')
        if not mummerplot_path or not gnuplot_path:
            return

        try:
            with open(mummerplot_path) as f:
                content = f.read()
        except (UnicodeDecodeError, PermissionError):
            return  # not the perl script we know how to fix

        if "$GNUPLOT_EXE = 'false'" not in content:
            return

        try:
            with open(mummerplot_path + '.tmp', 'w') as out_fh:
                for line in content.splitlines(keepends=True):
                    if "$GNUPLOT_EXE = 'false'" in line:
                        line = f"my $GNUPLOT_EXE = '{gnuplot_path}';\n"
                    out_fh.write(line)
            shutil.move(mummerplot_path + '.tmp', mummerplot_path)
            os.chmod(mummerplot_path, 0o0755)
        except PermissionError:
            sys.stderr.write('Could not fix the gnuplot path in mummerplot (no write permission): '
                             'the dot plots may be skipped.\n')
