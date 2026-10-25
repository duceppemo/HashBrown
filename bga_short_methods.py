import os
import sys
import shutil
import glob
import pathlib
import textwrap
import subprocess
from concurrent import futures
from psutil import virtual_memory
from multiprocessing import cpu_count


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz']

    @staticmethod
    def check_input(my_input):
        # Check if exists and is a folder
        if not (os.path.exists(my_input) and os.path.isdir(my_input)):
            raise Exception('Please select an existing folder as input.')

        # List content of folder
        file_list = os.listdir(my_input)

        # Check if folder is not empty and all files have the accepted extensions
        if not any([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure some files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_fastq_files(in_folder):
        sample_dict = dict()
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].split('_')[0]
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]

                    if sample not in sample_dict:
                        sample_dict[sample] = list()
                    if '_R1' in filename:
                        sample_dict[sample].insert(0, file_path)
                    elif '_R2' in filename:
                        sample_dict[sample].insert(1, file_path)

        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob.glob(folder + '/*' + extension)

    @staticmethod
    def format_fasta(input_fasta, output_fasta):
        fasta_dict = dict()

        # Parse fasta file into dictionary
        with open(input_fasta, 'r') as in_handle:
            header = ''
            seq = list()  # Use list to store sequence information
            for line in in_handle:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seq:
                        fasta_dict[header] = ''.join(seq)  # Store in dictionary
                        seq = list()  # Empty seq
                    header = line
                else:
                    seq.append(line)
                # For the last entry
                fasta_dict[header] = ''.join(seq)

        # Write to file with a width of 80 character max for the sequence
        # Any additional sequences will be written of the next line
        with open(output_fasta, 'w') as out_handle:
            for title, seq in fasta_dict.items():
                out_handle.write('>{}\n{}\n'.format(title, '\n'.join(textwrap.wrap(seq, 80, break_long_words=True))))

    @staticmethod
    def read_qc_falco(sample, path_list, output_folder, cpu):
        for read_file in path_list:
            cmd = ['falco',
                   '--outdir', output_folder,
                   '--threads', str(cpu),
                   '-report-filename', output_folder + os.path.basename(read_file).split('.')[0] + '.html',
                   '-skip-data',
                   '-skip-summary',
                   read_file]

            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def read_qc_falco_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, output_folder, int(cpu / parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.read_qc_falco(*x), args):
                pass

    @staticmethod
    def trim_fastp(sample, path_list, output_folder, report_folder, cpu):
        # I/O
        in_r1 = path_list[0]
        out_r1 = output_folder + sample + '_R1.fastq.gz'

        cmd = ['fastp',
               '--in1', in_r1,
               '--out1', out_r1,
               '--cut_right',
               '--cut_right_mean_quality', str(10),
               '--length_required', str(64),
               '--html', report_folder + sample + '.html',
               '--thread', str(cpu)]
        if len(path_list) == 2:  # paired-end data
            in_r2 = path_list[1]
            out_r2 = output_folder + sample + '_R2.fastq.gz'

            cmd += ['--in2', in_r2,
                    '--out2', out_r2,
                    '--detect_adapter_for_pe']

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def trim_fastp_parallel(sample_dict, output_folder, report_folder, cpu, parallel):
        Methods.make_folder(output_folder)
        Methods.make_folder(report_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, output_folder, report_folder, int(cpu / parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.trim_fastp(*x), args):
                pass

    @staticmethod
    def assemble_skesa(sample, path_list, output_folder, cpu, mem):
        # I/O
        r1 = path_list[0]
        assembly = output_folder + sample + '.fasta'

        cmd = ['skesa',
               '--cores', str(cpu),
               '--memory', str(mem),
               '--contigs_out', assembly,
               '--reads']
        if len(path_list) == 1:
            cmd += [r1]
        else:  # elif len(path_list) == 2:
            r2 = path_list[1]
            cmd += ['{},{}'.format(r1, r2)]

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def assemble_skesa_parallel(sample_dict, output_folder, cpu, mem, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, output_folder, int(cpu / parallel), mem)
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_skesa(*x), args):
                pass

    @staticmethod
    def assemble_spades(sample, path_list, output_folder, cpu, mem):
        # I/O
        sub_folder = output_folder + sample + '/'
        Methods.make_folder(sub_folder)
        r1 = path_list[0]

        cmd = ['spades.py',
               '--threads', str(cpu),
               '--memory', str(mem),
               '-o', sub_folder]
        if len(path_list) == 1:
            cmd += ['--s', str(1), r1]
        else:
            r2 = path_list[1]
            cmd += ['--pe-1', str(1), r1,
                    '--pe-2', str(1), r2]

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Move and rename assembly and assembly graph files
        try:  # In case there is no assembly
            shutil.move(sub_folder + 'scaffolds.fasta', output_folder + sample + '.fasta')
            shutil.move(sub_folder + 'assembly_graph_with_scaffolds.gfa', output_folder + sample + '.gfa')
        except FileNotFoundError:
            pass

        # Delete subfolder
        shutil.rmtree(sub_folder, ignore_errors=False, onerror=None)  # delete all

    @staticmethod
    def assemble_spades_parallel(sample_dict, output_folder, cpu, mem, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, output_folder, int(cpu / parallel), int(mem / parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_spades(*x), args):
                pass

    @staticmethod
    def gfa_connector(sample, path_list, assembly, output_folder, cpu, gfa_con):
        r1 = path_list[0]

        cmd = [gfa_con,
               '--cores', str(cpu),
               '--contigs', assembly,
               '--gfa', output_folder + sample + '.gfa',
               '--reads']
        if len(path_list) == 1:
            cmd += [r1]
        else:  # elif len(path_list) == 2:
            r2 = path_list[1]
            cmd += ['{},{}'.format(r1, r2)]

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def gfa_connector_parallel(assembly_list, sample_dict, output_folder, cpu, parallel, gfa_con):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, [assembly for assembly in assembly_list if sample in assembly][0],
                     output_folder, int(cpu / parallel), gfa_con) for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.gfa_connector(*x), args):
                pass

    @staticmethod
    def assembly_graph(gfa, output_folder):
        sample = os.path.basename(gfa).split('.')[0]
        cmd_bandage = ['Bandage', 'image', gfa, output_folder + sample + '_graph.png']
        subprocess.run(cmd_bandage, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def assembly_graph_parallel(gfa_list, output_folder, cpu):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(cpu)) as executor:
            args = ((gfa, output_folder) for gfa in gfa_list)
            for results in executor.map(lambda x: Methods.assembly_graph(*x), args):
                pass

    @staticmethod
    def map_minimap2(ref, read_folder, cpu, output_folder):
        # I/O
        sample = os.path.basename(ref).split('.')[0]
        r1 = read_folder + sample + '_R1.fastq.gz'
        r2 = read_folder + sample + '_R2.fastq.gz'
        output_bam = output_folder + sample + '.bam'

        # cmd_bwa_mem = ['bwa', 'mem', '-t', str(cpu), genome, r1, r2]
        minimap2_cmd = ['minimap2', '-a', '-x', 'sr', '-t', str(cpu), ref, r1]
        if os.path.exists(r2):
            minimap2_cmd += [r2]
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-F', '4', '-h', '-']
        samtools_fixmate_cmd = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '-']
        samtools_markdup_cmd = ['samtools', 'markdup', '-r', '-@', str(cpu), '-', output_bam]
        samtools_index_cmd = ['samtools', 'index', output_bam]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        if os.path.exists(r2):
            p3 = subprocess.Popen(samtools_fixmate_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2.stdout.close()
            p4 = subprocess.Popen(samtools_sort_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p3.stdout.close()
            p5 = subprocess.Popen(samtools_markdup_cmd, stdin=p4.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p4.stdout.close()
            p5.communicate()
        else:
            p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2.stdout.close()
            p4 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p3.stdout.close()
            p4.communicate()

        # Index bam file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_minimap2_parallel(assembly_list, read_folder, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((ref, read_folder, int(cpu / parallel), output_folder)
                    for ref in assembly_list)
            for results in executor.map(lambda x: Methods.map_minimap2(*x), args):
                pass

    @staticmethod
    def run_qualimap(bam, output_folder, cpu, mem):
        # I/O
        sample = os.path.basename(bam).split('.')[0]
        sub_folder = output_folder + sample + '/'
        Methods.make_folder(sub_folder)

        cmd = ['qualimap', 'bamqc',
               # '--paint-chromosome-limits',  # not so great with short read assemblies which result in many contigs
               '-bam', bam,
               '--java-mem-size={}G'.format(mem),
               '-nt', str(cpu),
               '-outdir', sub_folder,
               '-outformat', 'HTML']
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Rename html report
        os.rename(sub_folder + 'qualimapReport.html', sub_folder + sample + '.html')

    @staticmethod
    def run_qualimap_parallel(bam_list, output_folder, cpu, mem, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((bam, output_folder, int(cpu / parallel), int(mem / parallel))
                    for bam in bam_list)
            for results in executor.map(lambda x: Methods.run_qualimap(*x), args):
                pass
        
        # Cleanup bam files
        ext = ['.bam', '.bai']
        for i in ext:
            file_list = glob.glob(output_folder + '*' + i)
            for j in file_list:
                if os.path.exists(j):
                    os.remove(j)

    @staticmethod
    def run_quast(assembly_list, output_folder, cpu):
        cmd = ['quast.py',
               '--output-dir', output_folder,
               '--threads', str(cpu),
               '--min-contig', str(1000),
               '--no-icarus',
               ]
        cmd += assembly_list  # All the genome files
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
