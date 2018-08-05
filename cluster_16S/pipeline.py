"""
Design principles:
  1. reading this script should not require a high level of Python experience
  2. pipeline steps should correspond directly to functions
  3. there should be a single function that represents the pipeline
  4. for development the pipeline should be 'restartable'
"""
import argparse
import glob
import gzip
import itertools
import logging
import os
import re
import shutil
import sys

from cluster_16S.pipeline_util import create_output_dir, get_forward_fastq_files, get_associated_reverse_fastq_fp, \
    gzip_files, ungzip_files, run_cmd, PipelineException
from cluster_16S.fasta_qual_to_fastq import fasta_qual_to_fastq


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()

    Pipeline(**args.__dict__).run(input_dir=args.input_dir)
    return 0


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', default='.',
                            help='path to the input directory')
    arg_parser.add_argument('-w', '--work-dir', default='.',
                            help='path to the output directory')

    arg_parser.add_argument('-c', '--core-count', default=1, type=int,
                            help='number of cores to use')

    arg_parser.add_argument('--forward-primer', default='ATTAGAWACCCVNGTAGTCC',
                            help='forward primer to be clipped')
    arg_parser.add_argument('--reverse-primer', default='TTACCGCGGCKGCTGGCAC',
                            help='reverse primer to be clipped')

    arg_parser.add_argument('--uchime-ref-db-fp', default='/16SrDNA/pr2/pr2_gb203_version_4.5.fasta',
                            help='database for vsearch --uchime_ref')

    arg_parser.add_argument('--cutadapt-min-length', type=int, default=-1,
                            help='min_length for cutadapt, filling this in indicates that there are primers/adadpters to be removed')

    arg_parser.add_argument('--pear-min-overlap', required=True, type=int,
                            help='-v/--min-overlap for pear')
    
    arg_parser.add_argument('--pear-max-assembly-length', required=True, type=int,
                            help='-m/--max-assembly-length for pear')

    arg_parser.add_argument('--pear-min-assembly-length', required=True, type=int,
                            help='-m/--min-assembly-length for pear')

    arg_parser.add_argument('--vsearch-filter-maxee', required=True, type=int,
                            help='fastq_maxee for vsearch')
    
    arg_parser.add_argument('--vsearch-filter-trunclen', required=True, type=int,
                            help='fastq_trunclen for vsearch')

    arg_parser.add_argument('--vsearch-derep-minuniquesize', required=True, type=int,
                            help='minimum unique size for vsearch -derep_fulllength')

    args = arg_parser.parse_args()
    return args


class Pipeline:
    def __init__(self,
            work_dir,
            core_count,
            cutadapt_min_length,
            forward_primer, reverse_primer,
            pear_min_overlap, pear_max_assembly_length, pear_min_assembly_length,
            vsearch_filter_maxee, vsearch_filter_trunclen,
            vsearch_derep_minuniquesize,
            uchime_ref_db_fp,
            **kwargs  # allows some command line arguments to be ignored
    ):

        self.work_dir = work_dir
        self.core_count = core_count

        self.cutadapt_min_length = cutadapt_min_length
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

        self.pear_min_overlap = pear_min_overlap
        self.pear_max_assembly_length = pear_max_assembly_length
        self.pear_min_assembly_length = pear_min_assembly_length

        self.vsearch_filter_maxee = vsearch_filter_maxee
        self.vsearch_filter_trunclen = vsearch_filter_trunclen

        self.vsearch_derep_minuniquesize = vsearch_derep_minuniquesize

        self.uchime_ref_db_fp = uchime_ref_db_fp

        self.cutadapt_executable_fp = os.environ.get('CUTADAPT', default='cutadapt')
        self.pear_executable_fp = os.environ.get('PEAR', default='pear')
        self.usearch_executable_fp = os.environ.get('USEARCH', default='usearch')
        self.vsearch_executable_fp = os.environ.get('VSEARCH', default='vsearch')

    def run(self, input_dir):
        output_dir_list = list()
        output_dir_list.append(self.step_01_copy_and_compress(input_dir=input_dir))
        if self.cutadapt_min_length != -1:
            output_dir_list.append(self.step_01.5_remove_primers(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_02_merge_forward_reverse_reads_with_pear(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_qc_reads_with_vsearch(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_04_combine_runs(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_05_dereplicate_sort_remove_low_abundance_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_06_cluster_97_percent(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_07_reference_based_chimera_detection(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_08_create_otu_table(input_dir=output_dir_list[-1]))

        return output_dir_list

    def initialize_step(self):
        function_name = sys._getframe(1).f_code.co_name
        log = logging.getLogger(name=function_name)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        return log, output_dir

    def complete_step(self, log, output_dir):
        output_dir_list = sorted(os.listdir(output_dir))
        if len(output_dir_list) == 0:
            raise PipelineException('ERROR: no output files in directory "{}"'.format(output_dir))
        else:
            log.info('output files:\n\t%s', '\n\t'.join(os.listdir(output_dir)))
            # apply FastQC to all .fastq files
            fastq_file_list = [
                os.path.join(output_dir, output_file)
                for output_file
                in output_dir_list
                if re.search(pattern=r'\.fastq(\.gz)?$', string=output_file)
            ]

            if len(fastq_file_list) == 0:
                log.info('no .fastq files found in "{}"'.format(output_dir))
            else:
                fastqc_output_dir = os.path.join(output_dir, 'fastqc_results')
                os.makedirs(fastqc_output_dir, exist_ok=True)
                run_cmd(
                    [
                        'fastqc',
                        '--threads', str(self.core_count),
                        '--outdir', fastqc_output_dir,
                        *fastq_file_list
                    ],
                    log_file=os.path.join(fastqc_output_dir, 'log')
                )

    def step_01_copy_and_compress(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.debug('output_dir: %s', output_dir)
            # Check for fasta and qual files
            fasta_file_glob = os.path.join(input_dir, '*fasta*')
            qual_file_glob = os.path.join(input_dir, '*.qual*')
            fasta_files = sorted(glob.glob(fasta_file_glob))
            qual_files = sorted(glob.glob(qual_file_glob))
            fastas_no_slashes = [fasta.replace('\\', ' ').replace('/', ' ').split()[-1] for fasta in fasta_files]
            quals_no_slashes = [qual.replace('\\', ' ').replace('/', ' ').split()[-1] for qual in qual_files]
            fastas = [fasta.split(".")[0] for fasta in fastas_no_slashes]
            quals = [qual.split('.')[0] for qual in quals_no_slashes]
            for fasta in fastas:
                for qual in quals:
                    if fasta == qual:
                        tmpfasta = os.path.join(input_dir, fasta + ".fasta")
                        tmpqual = os.path.join(input_dir, qual + ".qual")
                        tmpfastq = os.path.join(input_dir, fasta + ".fastq")
                        log.info("Creating %s", tmpfastq)
                        fasta_qual_to_fastq(tmpfasta, tmpqual, tmpfastq)
                        log.info('%s created', tmpfastq)
                        quals.remove(qual)
                        break
            input_file_glob = os.path.join(input_dir, '*.fastq*')
            log.debug('input_file_glob: %s', input_file_glob)
            input_fp_list = sorted(glob.glob(input_file_glob))
            log.info('input files: %s', input_fp_list)

            if len(input_fp_list) == 0:
                raise PipelineException('found no fastq files in directory "{}"'.format(input_dir))

            for input_fp in input_fp_list:
                destination_fp = os.path.join(output_dir, os.path.basename(input_fp))
                if input_fp.endswith('.gz'):
                    with open(input_fp, 'rb') as f, open(destination_fp, 'wb') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)
                else:
                    destination_fp = destination_fp + '.gz'
                    with open(input_fp, 'rt') as f, gzip.open(destination_fp, 'wt') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)

        self.complete_step(log, output_dir)
        return output_dir
    
    
    def step_01.5_remove_primers(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('using cutadapt "%s"', self.cutadapt_executable_fp)

            for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                log.info('removing forward primers from file "%s"', forward_fastq_fp)
                forward_fastq_basename = os.path.basename(forward_fastq_fp)

                reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)
                log.info('removing reverse primers from file "%s"', reverse_fastq_fp)

                trimmed_forward_fastq_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=forward_fastq_basename,
                        pattern='_([0R])1',
                        repl=lambda m: '_trimmed_{}1'.format(m.group(1))))
                trimmed_reverse_fastq_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=forward_fastq_basename,
                        pattern='_([0R])1',
                        repl=lambda m: '_trimmed_{}2'.format(m.group(1))))

                run_cmd([
                        self.cutadapt_executable_fp,
                        '-a', self.forward_primer,
                        '-A', self.reverse_primer,
                        '-o', trimmed_forward_fastq_fp,
                        '-p', trimmed_reverse_fastq_fp,
                        '-m', str(self.cutadapt_min_length),
                        forward_fastq_fp,
                        reverse_fastq_fp
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                self.complete_step(log, output_dir)
        return output_dir
        

    """
    def step_03_merge_forward_reverse_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('vsearch executable: "%s"', self.vsearch_executable_fp)

            for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_merged'.format(m.group(1))
                )
                joined_fastq_fp = os.path.join(output_dir, joined_fastq_basename)[:-3]
                log.info('writing joined paired-end reads to "%s"', joined_fastq_fp)

                notmerged_fwd_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_notmerged_fwd'.format(m.group(1))
                )
                notmerged_fwd_fastq_fp = os.path.join(output_dir, notmerged_fwd_fastq_basename)[:-3]

                notmerged_rev_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_notmerged_rev'.format(m.group(1))
                )
                notmerged_rev_fastq_fp = os.path.join(output_dir, notmerged_rev_fastq_basename)[:-3]

                run_cmd([
                        self.vsearch_executable_fp,
                        '--fastq_mergepairs', forward_fastq_fp,
                        '--reverse', reverse_fastq_fp,
                        '--fastqout', joined_fastq_fp,
                        '--fastqout_notmerged_fwd', notmerged_fwd_fastq_fp,
                        '--fastqout_notmerged_rev', notmerged_rev_fastq_fp,
                        '--fastq_minovlen', str(self.pear_min_overlap),
                        '--fastq_maxlen', str(self.pear_max_assembly_length),
                        '--fastq_minlen', str(self.pear_min_assembly_length),
                        '--threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                gzip_files(glob.glob(os.path.join(output_dir, '*.fastq')))

        self.complete_step(log, output_dir)
        return output_dir
    """    


    def step_02_merge_forward_reverse_reads_with_pear(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            for compressed_forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                compressed_reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=compressed_forward_fastq_fp)

                forward_fastq_fp, reverse_fastq_fp = ungzip_files(
                    compressed_forward_fastq_fp,
                    compressed_reverse_fastq_fp,
                    target_dir=output_dir
                )

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_merged'.format(m.group(1)))[:-6]

                joined_fastq_fp_prefix = os.path.join(output_dir, joined_fastq_basename)
                log.info('joining paired ends from "%s" and "%s"', forward_fastq_fp, reverse_fastq_fp)
                log.info('writing joined paired-end reads to "%s"', joined_fastq_fp_prefix)
                run_cmd([
                        self.pear_executable_fp,
                        '-f', forward_fastq_fp,
                        '-r', reverse_fastq_fp,
                        '-o', joined_fastq_fp_prefix,
                        '--min-overlap', str(self.pear_min_overlap),
                        '--max-assembly-length', str(self.pear_max_assembly_length),
                        '--min-assembly-length', str(self.pear_min_assembly_length),
                        '-j', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                # delete the uncompressed input files
                os.remove(forward_fastq_fp)
                os.remove(reverse_fastq_fp)

                gzip_files(glob.glob(joined_fastq_fp_prefix + '.*.fastq'))

        self.complete_step(log, output_dir)
        return output_dir

    def step_03_qc_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_files_glob = os.path.join(input_dir, '*.assembled.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            for assembled_fastq_fp in glob.glob(input_files_glob):
                input_file_basename = os.path.basename(assembled_fastq_fp)
                output_file_basename = re.sub(
                    string=input_file_basename,
                    pattern='\.fastq\.gz',
                    repl='.ee{}trunc{}.fastq.gz'.format(self.vsearch_filter_maxee, self.vsearch_filter_trunclen)[:-3]
                )
                output_fastq_fp = os.path.join(output_dir, output_file_basename)

                log.info('vsearch executable: "%s"', self.vsearch_executable_fp)
                log.info('filtering "%s"', assembled_fastq_fp)
                run_cmd([
                        self.vsearch_executable_fp,
                        '-fastq_filter', assembled_fastq_fp,
                        '-fastqout', output_fastq_fp,
                        '-fastq_maxee', str(self.vsearch_filter_maxee),
                        '-fastq_trunclen', str(self.vsearch_filter_trunclen),
                        '-threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

            gzip_files(glob.glob(os.path.join(output_dir, '*.assembled.*.fastq')))

        self.complete_step(log, output_dir)
        return output_dir

    def step_04_combine_runs(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.assembled.*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))
            log.info('combining files:\n\t%s', '\n\t'.join(input_fp_list))

            output_file_name = get_combined_file_name(input_fp_list=input_fp_list)

            log.info('combined file: "%s"', output_file_name)
            output_fp = os.path.join(output_dir, output_file_name)
            with gzip.open(output_fp, 'wt') as output_file:
                for input_fp in input_fp_list:
                    with gzip.open(input_fp, 'rt') as input_file:
                        shutil.copyfileobj(fsrc=input_file, fdst=output_file)

        self.complete_step(log, output_dir)
        return output_dir

    def step_05_dereplicate_sort_remove_low_abundance_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.assembled.*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))

            for input_fp in input_fp_list:
                output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz$',
                        repl='.derepmin{}.fasta'.format(self.vsearch_derep_minuniquesize)))

                uc_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz$',
                        repl='.derepmin{}.txt'.format(self.vsearch_derep_minuniquesize)))

                run_cmd([
                        self.vsearch_executable_fp,
                        '-derep_fulllength', input_fp,
                        '-output', output_fp,
                        '-uc', uc_fp,
                        '-sizeout',
                        '-minuniquesize', str(self.vsearch_derep_minuniquesize),
                        '-threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

            gzip_files(glob.glob(os.path.join(output_dir, '*.fasta')))

        self.complete_step(log, output_dir)
        return output_dir

    def step_06_cluster_97_percent(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            # can usearch read gzipped files? no
            input_files_glob = os.path.join(input_dir, '*.fasta.gz')
            input_fp_list = sorted(glob.glob(input_files_glob))

            for compressed_input_fp in input_fp_list:

                input_fp, *_ = ungzip_files(compressed_input_fp, target_dir=output_dir)
                log.debug('input_fp: "%s"', input_fp)
                otu_output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.rad3.fasta'
                    )
                )

                uparse_output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.rad3.txt'
                    )
                )

                run_cmd([
                        self.usearch_executable_fp,
                        '-cluster_otus', input_fp,
                        '-otus', otu_output_fp,
                        '-relabel', 'OTU_',
                        # '-sizeout',
                        '-uparseout', uparse_output_fp
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                os.remove(input_fp)

        self.complete_step(log, output_dir)
        return output_dir

    def step_07_reference_based_chimera_detection(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_fps = glob.glob(os.path.join(input_dir, '*.fasta'))
            for input_fp in input_fps:
                uchimeout_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.uchime.txt'))

                notmatched_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.uchime.fasta'))

                log.info('starting chimera detection on file "%s"', input_fp)
                
                '''
                run_cmd([
                        self.usearch_executable_fp,
                        '-uchime2_ref', input_fp,
                        '-db', self.uchime_ref_db_fp,
                        '-uchimeout', uchimeout_fp,
                        '-mode', 'balanced',
                        '-strand', 'plus',
                        '-notmatched', notmatched_fp,
                        '-threads', str(self.core_count)
                    ],
                    log_file=os.path.join(output_dir, 'log')
                )
                '''

                run_cmd([
                        self.vsearch_executable_fp,
                        '-uchime_ref', input_fp,
                        '-db', self.uchime_ref_db_fp,
                        '-uchimeout', uchimeout_fp,
                        '-nonchimeras', notmatched_fp,
                        '-threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

        self.complete_step(log, output_dir)
        return output_dir

    def step_08_create_otu_table(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            otus_fp, *_ = glob.glob(os.path.join(input_dir, '*rad3.uchime.fasta'))
            input_fps = glob.glob(os.path.join(self.work_dir, 'step_03*', '*.assembled.fastq.gz'))
            for input_fp in input_fps:
                fasta_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz',
                        repl='.fasta'))

                log.info('convert fastq file\n\t%s\nto fasta file\n\t%s', input_fp, fasta_fp)
                run_cmd([
                    self.vsearch_executable_fp,
                    '--fastq_filter', input_fp,
                    '--fastaout', fasta_fp
                    ],
                    log_file=os.path.join(output_dir, 'log')
                )

                otu_table_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.assembled\.fastq\.gz$',
                        repl='.uchime.otutab.txt'
                    )
                )
                otu_table_biom_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.assembled\.fastq\.gz$',
                        repl='.uchime.otutab.json'
                    )
                )

                run_cmd([
                        self.vsearch_executable_fp,
                        '--usearch_global', fasta_fp,
                        '--db', otus_fp,
                        '--id', '0.97',
                        '--biomout', otu_table_biom_fp,
                        '--otutabout', otu_table_fp
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                # os.remove(fasta_fp)

        self.complete_step(log, output_dir)
        return output_dir

def get_combined_file_name(input_fp_list):
    if len(input_fp_list) == 0:
        raise PipelineException('get_combined_file_name called with empty input')

    def sorted_unique_elements(elements):
        return sorted(set(elements))

    return '_'.join(  # 'Mock_Run3_Run4_V4.fastq.gz'
        itertools.chain.from_iterable(  # ['Mock', 'Run3', 'Run4', 'V4.fastq.gz']
            map(  # [{'Mock'}, {'Run3', 'Run4'}, {'V4.fastq.gz'}]
                sorted_unique_elements,
                zip(  # [('Mock', 'Mock'), ('Run4', 'Run3'), ('V4.fastq.gz', 'V4.fastq.gz')]
                    *[  # [('Mock', 'Run3', 'V4.fastq.gz'), ('Mock', 'Run4', 'V4.fastq.gz')]
                        os.path.basename(fp).split('_')  # ['Mock', 'Run3', 'V4.fastq.gz']
                        for fp  # '/some/data/Mock_Run3_V4.fastq.gz'
                        in input_fp_list  # ['/input/data/Mock_Run3_V4.fastq.gz', '/input_data/Mock_Run4_V4.fastq.gz']
                    ]))))


if __name__ == '__main__':
    main()
