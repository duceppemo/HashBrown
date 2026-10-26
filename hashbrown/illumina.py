"""Short-read steps: trimming (fastp), assembly (SKESA or SPAdes) and polishing (Polypolish, pyPOLCA)."""

import glob
import io
import os
import shutil
import subprocess
from concurrent import futures

from hashbrown.methods import Methods


class IlluminaMethods:
    @staticmethod
    def trim_fastp(sample, sample_obj, output_folder, report_folder, cpu, flag, ion_torrent=False):
        """Trim a single-end file or a pair with fastp. Returns the trimmed file list.

        Ion Torrent reads keep the sliding-window quality trim (their quality decays toward the 3' end
        too) but skip adapter trimming: the Torrent Suite already removed the adapters, and fastp's
        adapter detection is tuned to Illumina chemistry.
        """
        reads = sample_obj.illumina.raw
        trimmed = list()

        if not reads:  # sample without short reads (hybrid run with missing pairs)
            return sample, trimmed

        if len(reads) == 2:
            trimmed = [os.path.join(output_folder, sample + '_R1.fastq.gz'),
                       os.path.join(output_folder, sample + '_R2.fastq.gz')]
        else:
            trimmed = [os.path.join(output_folder, sample + '.fastq.gz')]

        if not os.path.exists(flag):
            cmd = ['fastp',
                   '--in1', reads[0],
                   '--out1', trimmed[0],
                   '--cut_right',
                   '--cut_right_mean_quality', str(10),
                   '--length_required', str(64),
                   '--html', os.path.join(report_folder, sample + '.html'),
                   '--json', os.path.join(report_folder, sample + '.json'),
                   '--thread', str(cpu)]
            if len(reads) == 2:
                cmd += ['--in2', reads[1],
                        '--out2', trimmed[1]]
                if not ion_torrent:
                    cmd += ['--detect_adapter_for_pe']
            if ion_torrent:
                cmd += ['--disable_adapter_trimming']

            print(f'\t{sample}')
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        if not all(os.path.exists(f) for f in trimmed):
            trimmed = list()

        return sample, trimmed

    @staticmethod
    def trim_fastp_parallel(sample_dict, output_folder, report_folder, cpu, parallel, flag,
                            ion_torrent=False):
        Methods.make_folder(output_folder)
        Methods.make_folder(report_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, report_folder, max(1, int(cpu / parallel)), flag,
                     ion_torrent)
                    for sample, sample_obj in sample_dict.items())
            for sample, trimmed in executor.map(lambda x: IlluminaMethods.trim_fastp(*x), args):
                sample_dict[sample].illumina.trimmed = trimmed

        if os.path.exists(flag):
            print('\tSkipping trimming short reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def assemble_skesa(sample, sample_obj, output_folder, gfa_folder, cpu, mem, flag):
        reads = sample_obj.best_illumina()
        assembly = os.path.join(output_folder, sample + '.fasta')

        if not os.path.exists(flag) and reads:
            print(f'\t{sample}')
            cmd = ['skesa',
                   '--cores', str(cpu),
                   '--memory', str(mem),
                   '--contigs_out', assembly,
                   '--reads', ','.join(reads)]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # SKESA has no graph output: gfa_connector (from the SKESA authors) rebuilds one
            if os.path.exists(assembly) and shutil.which('gfa_connector'):
                Methods.make_folder(gfa_folder)
                gfa_file = os.path.join(gfa_folder, sample + '.gfa')
                cmd_gfa = ['gfa_connector',
                           '--cores', str(cpu),
                           '--contigs', assembly,
                           '--gfa', gfa_file,
                           '--reads', ','.join(reads)]
                subprocess.run(cmd_gfa, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                if os.path.exists(gfa_file):
                    Methods.make_assembly_graph_image(gfa_file, os.path.join(gfa_folder, sample + '.png'))

        if not os.path.exists(assembly):
            print(f'\tAssembly failed for {sample}')
            assembly = ''

        return sample, assembly

    @staticmethod
    def assemble_spades(sample, sample_obj, output_folder, gfa_folder, ion_torrent, cpu, mem, flag):
        reads = sample_obj.best_illumina()
        assembly = os.path.join(output_folder, sample + '.fasta')
        output_subfolder = os.path.join(output_folder, sample)

        if not os.path.exists(flag) and reads:
            print(f'\t{sample}')
            cmd = ['spades.py',
                   '--threads', str(cpu),
                   '--memory', str(mem),
                   '-o', output_subfolder]
            if ion_torrent:
                cmd += ['--iontorrent']
            if len(reads) == 2:
                cmd += ['-1', reads[0], '-2', reads[1]]
            else:
                cmd += ['-s', reads[0]]

            Methods.make_folder(output_subfolder)
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # Keep the scaffolds and the assembly graph, drop the working folder
            if os.path.exists(os.path.join(output_subfolder, 'scaffolds.fasta')):
                shutil.move(os.path.join(output_subfolder, 'scaffolds.fasta'), assembly)

                gfa = os.path.join(output_subfolder, 'assembly_graph_with_scaffolds.gfa')
                if os.path.exists(gfa):
                    Methods.make_folder(gfa_folder)
                    gfa_file = os.path.join(gfa_folder, sample + '.gfa')
                    shutil.move(gfa, gfa_file)
                    Methods.make_assembly_graph_image(gfa_file, os.path.join(gfa_folder, sample + '.png'))
            shutil.rmtree(output_subfolder, ignore_errors=True)

        if not os.path.exists(assembly):
            print(f'\tAssembly failed for {sample}')
            assembly = ''

        return sample, assembly

    @staticmethod
    def assemble_short_parallel(sample_dict, assembler, output_folder, gfa_folder, ion_torrent, cpu, mem,
                                parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            if assembler == 'skesa':
                args = ((sample, sample_obj, output_folder, gfa_folder,
                         max(1, int(cpu / parallel)), mem, flag)
                        for sample, sample_obj in sample_dict.items())
                runner = IlluminaMethods.assemble_skesa
            else:  # spades
                args = ((sample, sample_obj, output_folder, gfa_folder, ion_torrent,
                         max(1, int(cpu / parallel)), max(1, int(mem / parallel)), flag)
                        for sample, sample_obj in sample_dict.items())
                runner = IlluminaMethods.assemble_spades
            for sample, assembly in executor.map(lambda x: runner(*x), args):
                sample_dict[sample].assembly.raw = assembly

        if os.path.exists(flag):
            print('\tSkipping assembling short reads. Already done.')
        else:
            Methods.flag_done(flag)

    @staticmethod
    def fix_fasta(input_fasta, output_fasta):
        """Normalize line length and uppercase the sequences (the polishers write lowercase edits)."""
        cmd = ['seqtk', 'seq', input_fasta]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        with open(output_fasta, 'w') as f:
            for line in io.TextIOWrapper(p.stdout, encoding='utf-8'):
                if line.startswith('>'):
                    f.write(line)
                else:
                    f.write(line.upper())

    @staticmethod
    def fix_fasta_header(genome):
        """Rename the contigs 'contig_1', 'contig_2', ... dropping the polishing suffixes."""
        tmp_file = genome + '.tmp'
        os.rename(genome, tmp_file)

        with open(genome, 'w') as out_fh, open(tmp_file) as in_fh:
            i = 0
            for line in in_fh:
                if line.startswith('>'):
                    i += 1
                    out_fh.write(f'>contig_{i}\n')
                else:
                    out_fh.write(line)

        os.remove(tmp_file)

    @staticmethod
    def run_polypolish(genome, r1, r2, polished_folder, sample, cpu):
        # I/O
        sam_r1 = os.path.join(polished_folder, sample + '_R1.sam')
        sam_r2 = os.path.join(polished_folder, sample + '_R2.sam')

        # Index the genome, then map each mate separately keeping all alignments (-a), as Polypolish wants
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        for sam, reads in ((sam_r1, r1), (sam_r2, r2)):
            p = subprocess.Popen(['bwa', 'mem', '-t', str(cpu), '-a', genome, reads],
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            with open(sam, 'wb') as f:
                f.write(p.communicate()[0])

        cmd_pp = ['polypolish', 'polish', genome, sam_r1, sam_r2]
        p = subprocess.Popen(cmd_pp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        fasta_pp = os.path.join(polished_folder, sample + '.polypolish.fasta')
        with open(fasta_pp, 'wb') as f:
            f.write(stdout)
        with open(os.path.join(polished_folder, sample + '.polypolish.log'), 'wb') as f:
            f.write(stderr)

        # Fix fasta
        fixed_fasta = os.path.join(polished_folder, sample + '.fasta')
        IlluminaMethods.fix_fasta(fasta_pp, fixed_fasta)

        # Cleanup, including the index files bwa wrote next to the input genome
        for ext in ('.amb', '.ann', '.bwt', '.pac', '.sa', '.sam'):
            for leftover in glob.glob(os.path.join(polished_folder, sample + '*' + ext)) \
                    + glob.glob(genome + ext):
                os.remove(leftover)
        os.remove(fasta_pp)

        return fixed_fasta

    @staticmethod
    def run_pypolca(genome, r1, r2, polished_folder, sample, cpu):
        # https://github.com/gbouras13/pypolca
        # The polished fasta is <prefix>_corrected.fasta in the output folder; the report is <prefix>.report
        tmp_folder = os.path.join(polished_folder, sample)
        fasta_polca = os.path.join(tmp_folder, sample + '_corrected.fasta')
        Methods.make_folder(tmp_folder)

        cmd = ['pypolca', 'run',
               '--assembly', genome,
               '--reads1', r1,
               '--reads2', r2,
               '--output', tmp_folder,
               '--prefix', sample,
               '--careful',
               '--force',
               '--threads', str(cpu)]

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, _ = p.communicate()
        with open(os.path.join(polished_folder, sample + '.pypolca.log'), 'wb') as f:
            f.write(stdout)

        # Fix fasta
        fixed_fasta = os.path.join(polished_folder, sample + '.fasta')
        if os.path.exists(fasta_polca):
            IlluminaMethods.fix_fasta(fasta_polca, fixed_fasta)

        # Keep the report and the corrections, drop the working folder
        for kept in (sample + '.report', sample + '.vcf'):  # the vcf is sometimes absent
            src = os.path.join(tmp_folder, kept)
            dst = os.path.join(polished_folder, kept)
            if os.path.exists(src):
                if os.path.exists(dst):
                    os.remove(dst)
                shutil.move(src, dst)
        shutil.rmtree(tmp_folder, ignore_errors=True)

        return fixed_fasta

    @staticmethod
    def polish(sample, sample_obj, output_folder, polishers, cpu, flag):
        """Run the chosen short-read polishers, in order, on the best long-read assembly.

        polishers is a sequence of 'polypolish' and/or 'pypolca'. Each polisher takes over the
        previous one's output. The default, Polypolish then pyPOLCA --careful, is the combination
        that introduced no errors in low-coverage regions in our tests (see the wiki).
        """
        reads = sample_obj.best_illumina()
        genome = sample_obj.assembly.medaka or sample_obj.assembly.raw
        polished_assembly = os.path.join(output_folder, sample + '.fasta')

        if not os.path.exists(flag):
            if not genome or not os.path.exists(genome):
                print(f'\tNo assembly for {sample}. Skipping short read polishing.')
            elif len(reads) != 2:
                print(f'\tPaired-end data missing for {sample}. Skipping short read polishing.')
            else:
                print(f'\t{sample}')
                Methods.make_folder(output_folder)

                runners = {'polypolish': IlluminaMethods.run_polypolish,
                           'pypolca': IlluminaMethods.run_pypolca}
                for polisher in polishers:
                    genome = runners[polisher](genome, reads[0], reads[1], output_folder, sample, cpu)

                if os.path.exists(polished_assembly):
                    IlluminaMethods.fix_fasta_header(polished_assembly)

        if not os.path.exists(polished_assembly):
            polished_assembly = ''

        return sample, polished_assembly

    @staticmethod
    def polish_parallel(sample_dict, output_folder, polishers, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, polishers, max(1, int(cpu / parallel)), flag)
                    for sample, sample_obj in sample_dict.items())
            for sample, polished in executor.map(lambda x: IlluminaMethods.polish(*x), args):
                sample_dict[sample].assembly.polypolish = polished

        if os.path.exists(flag):
            print('\tSkipping short read polishing. Already done.')
        else:
            Methods.flag_done(flag)
