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

#TODO CHANGE BACK TO cluster_16S.
from cluster_16S.pipeline_util import create_output_dir, get_forward_fastq_files, get_associated_reverse_fastq_fp, \
    gzip_files, ungzip_files, run_cmd, PipelineException
from cluster_16S.fasta_qual_to_fastq import fasta_qual_to_fastq


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    args = check_args(args, **args.__dict__)
    Pipeline(**args.__dict__).run(input_dir=args.input_dir)
    return 0


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', required=True,
                            help='path to the input directory')

    arg_parser.add_argument('-w', '--work-dir', required=False, default="",
                            help='path to the output directory')

    arg_parser.add_argument('-c', '--core-count', default=1, type=int,
                            help='number of cores to use')

    arg_parser.add_argument('-m', '--multiple-runs', action='store_true', default=False,
                            help='indicates whether samples are split across multiple runs. If so, files must be named \"SAMPLENAME_run1.fastq\", \"SAMPLENAME_run2.fastq\", ...')

    arg_parser.add_argument('-p', '--paired-ends', action='store_true', default=False,
                            help='flag that indicates whether your data is in paired ends. If so, files must be named \"SAMPLENAME_R1.fastq\", \"SAMPLENAME_R1.fastq\"')

    arg_parser.add_argument('--steps', default=6, type=int,
                            help='Indicates that pipeline should only run up to this step. Optional steps (such as remove_primers or combine_runs) will not count against this step count')

    arg_parser.add_argument('--debug', default=False, action='store_true',
                            help='Flag to start debugging output. Errors and warning will still be written to stderr if this flag is not used')

    arg_parser.add_argument('--forward-primer', default='ATTAGAWACCCVNGTAGTCC',
                            help='Forward primer to be clipped by cutadapt')

    arg_parser.add_argument('--reverse-primer', default='TTACCGCGGCKGCTGGCAC',
                            help='Reverse primer to be clipped by cutadapt')

    arg_parser.add_argument('--uchime-ref-db-fp', default='/app/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz',
                            help='Database for vsearch --uchime_ref (if using singularity image, will use built-in SILVA database if no arg given)')

    arg_parser.add_argument('--cutadapt-min-length', type=int, default=-1,
                            help='Min_length for cutadapt, filling this in indicates that there are primers/adadpters to be removed and cutadapt will be used')

    arg_parser.add_argument('--cutadapt-adapter-file-forward', default="",
                            help='path to the file containing all forward adapters in your samples')

    arg_parser.add_argument('--cutadapt-adapter-file-reverse', default="",
                            help='path to the file containing all reverse adapters in your samples')

    arg_parser.add_argument('--pear-min-overlap', required=False, type=int, default=10,
                            help='-v/--min-overlap for pear')
    
    arg_parser.add_argument('--pear-max-assembly-length', required=False, type=int, default=500,
                            help='-m/--max-assembly-length for pear')

    arg_parser.add_argument('--pear-min-assembly-length', required=False, type=int, default=200,
                            help='-m/--min-assembly-length for pear')

    arg_parser.add_argument('--vsearch-filter-maxee', required=True, type=int,
                            help='fastq_maxee for vsearch')
    
    arg_parser.add_argument('--vsearch-filter-trunclen', required=True, type=int,
                            help='fastq_trunclen for vsearch')

    arg_parser.add_argument('--vsearch-derep-minuniquesize', required=True, type=int,
                            help='Minimum unique size for vsearch -derep_fulllength')

    arg_parser.add_argument('--combine-final-results', action='store_true', default=False,
                            help='Flag that indicates that all samples should be combined into 1 OTU table') 
    
    args = arg_parser.parse_args()
    return args

def check_args(args,
            input_dir,
            work_dir,
            core_count,
            multiple_runs,
            paired_ends,
            steps,
            debug,
            cutadapt_min_length,
            cutadapt_adapter_file_forward,
            cutadapt_adapter_file_reverse,
            forward_primer, reverse_primer,
            pear_min_overlap, pear_max_assembly_length, pear_min_assembly_length,
            vsearch_filter_maxee, vsearch_filter_trunclen,
            vsearch_derep_minuniquesize,
            uchime_ref_db_fp,
            combine_final_results,
            **kwargs  # allows some command line arguments to be ignored
    ):
    if not os.path.isdir(input_dir):
        print("{} is not a directory".format(input_dir), file=sys.stderr)
        exit()
    if len(os.listdir(input_dir)) == 0:
        print("{} is empty".format(input_dir), file=sys.stderr)
        exit()
    if len([f for f in os.listdir(input_dir) if os.path.isfile(f)]) % 2 == 1:
        print("Odd number of input files in {}. Please check to make sure there are no missing paired-end or run files".format(input_dir), file=sys.stderr)
        exit()
    if work_dir == "":
        args.work_dir = "{}/output".format(input_dir)
        if not os.path.isdir(args.work_dir):
            os.mkdir(args.work_dir)
    if core_count < 1:
        print("Bad number of cores")
        exit()
    if cutadapt_min_length < -1 or cutadapt_min_length == 0:
        print("Invalid value for --cutadapt-min-length", file=sys.stderr)
        exit()
    if cutadapt_adapter_file_forward is not "" and not os.path.isfile(cutadapt_adapter_file_forward):
        print("{} is not a file".format(cutadapt_adapter_file_forward), file=sys.stderr)
        exit()
    if cutadapt_adapter_file_reverse is not "" and not os.path.isfile(cutadapt_adapter_file_reverse):
        print("{} is not a file".format(cutadapt_adapter_file_reverse), file=sys.stderr)
        exit()
    if paired_ends is True and ((cutadapt_adapter_file_forward is "" and cutadapt_adapter_file_reverse is not "") or (cutadapt_adapter_file_forward is not "" and cutadapt_adapter_file_reverse is "")):
        print("Paired ends is selected and only one adapter file is passed. You most likely forgot to include the other adapter file", file=sys.stderr)
        exit()
    if pear_min_overlap <= 0:
        print("Invalid value for --pear-min-overlap", file=sys.stderr)
        exit()
    if pear_max_assembly_length <= 0:
        print("Invalid value for --pear-max-assembly-length", file=sys.stderr)
        exit()
    if pear_min_assembly_length <= 0:
        print("Invalid value for --pear-min-assembly-length", file=sys.stderr)
        exit()
    if pear_min_assembly_length > pear_max_assembly_length:
        print("--pear-min-assembly-length can't be greater than --pear-max-assembly-length", file=sys.stderr)
        exit()
    if vsearch_filter_maxee < 0:
        print("Invalid value for --vsearch-filter-maxee", file=sys.stderr)
        exit()
    if vsearch_filter_trunclen <= 0:
        print("Invalid value for --vesearch-filter-trunclen", file=sys.stderr)
        exit()
    if vsearch_derep_minuniquesize < 0:
        print("Invalid value for --vsearch-derep-minuniquesize", file=sys.stderr)
        exit()
    if not os.path.isfile(uchime_ref_db_fp):
        print("{} is not a valid file path".format(uchime_ref_db_fp), file=sys.stderr)
        exit()
    if steps < 1:
        print("Invalid value for --steps", file=sys.stderr)
        exit()
    return args
 

class Pipeline:
    def __init__(self,
            work_dir,
            core_count,
            multiple_runs,
            paired_ends,
            steps,
            debug,
            cutadapt_min_length,
            cutadapt_adapter_file_forward,
            cutadapt_adapter_file_reverse,
            forward_primer, reverse_primer,
            pear_min_overlap, pear_max_assembly_length, pear_min_assembly_length,
            vsearch_filter_maxee, vsearch_filter_trunclen,
            vsearch_derep_minuniquesize,
            uchime_ref_db_fp,
            combine_final_results,
            **kwargs  # allows some command line arguments to be ignored
    ):

        self.work_dir = work_dir
        self.core_count = core_count
        self.multiple_runs = multiple_runs
        self.paired_ends = paired_ends
        self.steps = steps
        self.debug = debug
        self.cutadapt_min_length = cutadapt_min_length
        self.cutadapt_adapter_file_forward = cutadapt_adapter_file_forward
        self.cutadapt_adapter_file_reverse = cutadapt_adapter_file_reverse
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.combine_final_results = combine_final_results

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
        step_counter = 0
        output_dir_list.append(self.step_01_copy_and_compress(input_dir=input_dir))
        step_counter += 1
        if self.steps == step_counter:
            return output_dir_list
        if self.cutadapt_min_length != -1:
            output_dir_list.append(self.step_01_1_remove_primers(input_dir=output_dir_list[-1]))
        if self.paired_ends is True:
            output_dir_list.append(self.step_01_2_merge_forward_reverse_reads_with_pear(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_02_qc_reads_with_vsearch(input_dir=output_dir_list[-1]))
        step_counter += 1
        if self.steps == step_counter:
            return output_dir_list
        if self.multiple_runs is True:
            output_dir_list.append(self.step_02_1_combine_runs(input_dir=output_dir_list[-1]))
        if self.combine_final_results is True:
            output_dir_list.append(self.step_02_2_combine_samples(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_dereplicate_sort_remove_low_abundance_reads(input_dir=output_dir_list[-1]))
        step_counter += 1
        if self.steps == step_counter:
            return output_dir_list
        output_dir_list.append(self.step_04_cluster_97_percent(input_dir=output_dir_list[-1]))
        step_counter += 1
        if self.steps == step_counter:
            return output_dir_list
        output_dir_list.append(self.step_05_reference_based_chimera_detection(input_dir=output_dir_list[-1]))
        step_counter += 1
        if self.steps == step_counter:
            return output_dir_list
        output_dir_list.append(self.step_06_create_otu_table(input_dir=output_dir_list[-1]))

        return output_dir_list

    def initialize_step(self):
        function_name = sys._getframe(1).f_code.co_name
        log = logging.getLogger(name=function_name)
        if self.debug is True:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.WARNING)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir, debug=self.debug)
        return log, output_dir

    def complete_step(self, log, output_dir):
        output_dir_list = sorted(os.listdir(output_dir))
        if len(output_dir_list) == 0:
            raise PipelineException('ERROR: no output files in directory "{}"'.format(output_dir))
        return


    def step_01_copy_and_compress(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.debug('output_dir: %s', output_dir)
            # Check for fasta and qual files
            fasta_files = []
            types = ['*.fasta*', '*.fa*']
            for t in types:
                fasta_files.extend(glob.glob(os.path.join(input_dir, t)))
            fasta_files = sorted(fasta_files)
            log.info('fasta files: "%s"', fasta_files)
            qual_file_glob = os.path.join(input_dir, '*.qual*')
            qual_files = sorted(glob.glob(qual_file_glob))
            log.info('qual files: "%s"', qual_files)
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
            input_fp_list = []
            types = ['*.fastq*', '*.fq*']
            for t in types:
                input_fp_list.extend(glob.glob(os.path.join(input_dir, t)))
            input_fp_list = sorted(glob.glob(input_file_glob))
            log.info('input files: %s', input_fp_list)

            if len(input_fp_list) == 0:
                raise PipelineException('found no fastq files in directory "{}"'.format(input_dir))

            for input_fp in input_fp_list:
                destination_name = re.sub(
                                        string=os.path.basename(input_fp),
                                        pattern='\.fq',
                                        repl='.fastq')
                destination_fp = os.path.join(output_dir, destination_name)
                if input_fp.endswith('.gz'):
                    with open(input_fp, 'rb') as f, open(destination_fp, 'wb') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)
                else:
                    destination_fp = destination_fp + '.gz'
                    with open(input_fp, 'rt') as f, gzip.open(destination_fp, 'wt') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)

        self.complete_step(log, output_dir)
        return output_dir
    
    
    def step_01_1_remove_primers(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('using cutadapt "%s"', self.cutadapt_executable_fp)
            if self.paired_ends is True:
                for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir, debug=self.debug):
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
                    if self.cutadapt_adapter_file_forward is "":
                        run_cmd([
                                self.cutadapt_executable_fp,
                                '-a', self.forward_primer,
                                '-A', self.reverse_primer,
                                '-o', trimmed_forward_fastq_fp,
                                '-p', trimmed_reverse_fastq_fp,
                                '-m', str(self.cutadapt_min_length),
                                '-j', str(self.core_count),
                                forward_fastq_fp,
                                reverse_fastq_fp
                            ],
                            log_file = os.path.join(output_dir, 'log'),
                            debug=self.debug
                        )
                    else:
                        log.info('using adapters from file "%s" and "%s"', self.cutadapt_adapter_file_forward, self.cutadapt_adapter_file_reverse)
                        run_cmd([
                            self.cutadapt_executable_fp,
                            '-a', 'file:{}'.format(self.cutadapt_adapter_file_forward),
                            '-A', 'file:{}'.format(self.cutadapt_adapter_file_reverse),
                            '-o', trimmed_forward_fastq_fp,
                            '-p', trimmed_reverse_fastq_fp,
                            '-m', str(self.cutadapt_min_length),
                            '-j', str(self.core_count),
                            forward_fastq_fp,
                            reverse_fastq_fp
                        ],
                            log_file=os.path.join(output_dir, 'log'),
                            debug=self.debug
                        )
            else:
                input_files_glob = os.path.join(input_dir, '*.fastq.gz')
                input_file_list = glob.glob(input_files_glob)
                for input_file in input_file_list:
                    log.info('removing forward primers from file "%s"', input_file)
                    input_basename = os.path.basename(input_file)
                    trimmed_fastq_fp = os.path.join(
                                                output_dir,
                                                re.sub(
                                                    string=input_basename,
                                                    pattern='\.fastq\.gz$',
                                                    repl='_trimmed.fastq.gz'))
                    if self.cutadapt_adapter_file_forward is not "":
                        run_cmd([
                            self.cutadapt_executable_fp,
                            '-a', 'file:{}'.format(self.cutadapt_adapter_file_forward),
                            '-o', trimmed_fastq_fp,
                            '-m', str(self.cutadapt_min_length),
                            '-j', str(self.core_count),
                            input_file
                        ],
                            log_file=os.path.join(output_dir, 'log'),
                            debug=self.debug
                        )
                    else:
                        run_cmd([
                            self.cutadapt_executable_fp,
                            '-a', self.forward_primer,
                            '-o', trimmed_fastq_fp,
                            '-m', str(self.cutadapt_min_length),
                            '-j', str(self.core_count),
                            input_file
                        ],
                            log_file=os.path.join(output_dir, 'log'),
                            debug=self.debug
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


    def step_01_2_merge_forward_reverse_reads_with_pear(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            for compressed_forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir, debug=self.debug):
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
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

                # delete the uncompressed input files
                os.remove(forward_fastq_fp)
                os.remove(reverse_fastq_fp)
                gzip_files(glob.glob(joined_fastq_fp_prefix + '.*.fastq'), debug=self.debug)
                with open(os.path.join(output_dir, 'log'), 'r') as logcheck:
                    num_assembled = 0
                    num_discarded = 0
                    forward_fp = ""
                    reverse_fp = ""
                    for l in logcheck:
                        if 'Forward reads file' in l:
                            forward_fp = l.split(' ')[-1]
                        elif 'Reverse reads file' in l:
                            reverse_fp = l.split(' ')[-1]
                        elif 'Assembled reads' in l and 'file' not in l:
                            num_assembled = int(l.split(' ')[3].replace(',', ''))
                        elif 'Discarded reads' in l and 'file' not in l:
                            num_discarded = int(l.split(' ')[3].replace(',', ''))
                            log.info("num_assembled = {}, num_discarded = {}".format(num_assembled, num_discarded))
                            if num_discarded > num_assembled:
                                log.warning("More sequences discarded than kept by PEAR for files '{}' and '{}'".format(forward_fp, reverse_fp))

        self.complete_step(log, output_dir)
        return output_dir

    def step_02_qc_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            if self.paired_ends is True:
                input_files_glob = os.path.join(input_dir, '*.assembled*.fastq.gz')
            else:
                input_files_glob = os.path.join(input_dir, '*.fastq.gz')
            input_file_list = glob.glob(input_files_glob)
            if len(input_file_list) == 0:
                raise PipelineException('found no .fastq.gz files in directory "{}"'.format(input_dir))
            log.info('input file glob: "%s"', input_files_glob)
            for assembled_fastq_fp in input_file_list:
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
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
            gzip_files(glob.glob(os.path.join(output_dir, '*.fastq')), debug=self.debug)
            with open(os.path.join(output_dir, 'log'), 'r') as logcheck:
                kept_num = 0
                discarded_num = 0
                for l in logcheck:
                    if 'sequences kept' in l:
                        l_arr = l.split(' ')
                        kept_num = int(l_arr[0])
                        discarded_num = int(l_arr[7])
                    if 'executing' in l:
                        l_arr = l.split(' ')
                        ran_fp = l[3]
                        log.info("kept_num = {}, discarded_num = {}".format(kept_num, discarded_num))
                        if kept_num == 0:
                            log.error("No sequences kept by vsearch qc for input file '{}'".format(ran_fp))
                        if discarded_num > kept_num:
                            log.warning("More sequences discarded than kept by vsearch qc for input file '{}'".format(ran_fp))
        self.complete_step(log, output_dir)
        return output_dir

    def step_02_1_combine_runs(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*run1*.*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            run_fp_list = sorted(glob.glob(input_files_glob))
            if len(run_fp_list) == 0:
                raise PipelineException('found no run1*.fastq.gz files in directory "{}"'.format(input_dir))
            for run in run_fp_list:
                sample_name = os.path.basename(run).split('_run1')[0]
                trailing_name = os.path.basename(run).split('_run1')[1]
                log.info('Sample name: "%s"', sample_name)
                run_files_glob = os.path.join(input_dir, '%s*.fastq.gz' % sample_name)
                run_files_list = sorted(glob.glob(run_files_glob))
                log.info('Sample run file list: "%s"', run_files_list)
                output_run_file = '{}_combined{}'.format(sample_name, trailing_name)
                output_fp = os.path.join(output_dir, output_run_file)
                log.info('combined file: "%s"', output_fp)
                with gzip.open(output_fp, 'wt') as output_file:
                    for run_fp in run_files_list:
                        with gzip.open(run_fp, 'rt') as input_file:
                            shutil.copyfileobj(fsrc=input_file, fdst=output_file)
                    
        self.complete_step(log, output_dir)
        return output_dir

    def step_02_2_combine_samples(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.fastq.gz')
            sample_fp_list = sorted(glob.glob(input_files_glob))
            log.info('file list: {}'.format(sample_fp_list))
            if len(sample_fp_list) == 0:
                raise PipelineException('found no *.fastq.gz files in directory "{}"'.format(input_dir))
            combined_fp = '{}.fastq.gz'.format(self.get_combined_file_name(output_dir))
            log.info('combining into "{}"'.format(combined_fp))
            with gzip.open(combined_fp, 'wt') as output_file:
                for sample_fp in sample_fp_list:
                    with gzip.open(sample_fp, 'rt') as input_file:
                        shutil.copyfileobj(fsrc=input_file, fdst=output_file)
            self.complete_step(log, output_dir)
            return output_dir

    def step_03_dereplicate_sort_remove_low_abundance_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))
            if len(input_fp_list) == 0:
                raise PipelineException('found no .fastq.gz files in directory "{}"'.format(input_dir))
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
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

            gzip_files(glob.glob(os.path.join(output_dir, '*.fasta')), debug=self.debug)

        self.complete_step(log, output_dir)
        return output_dir

    def step_04_cluster_97_percent(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            # can usearch read gzipped files? no
            input_files_glob = os.path.join(input_dir, '*.fasta.gz')
            input_fp_list = sorted(glob.glob(input_files_glob))
            if len(input_fp_list) == 0:
                raise PipelineException('found no fasta.gz files in directory "{}"'.format(input_dir))
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
                        #'-sizeout',
                        '-uparseout', uparse_output_fp
                    ],
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

                os.remove(input_fp)

        self.complete_step(log, output_dir)
        return output_dir

    def step_05_reference_based_chimera_detection(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_fps = glob.glob(os.path.join(input_dir, '*.fasta'))
            if len(input_fps) == 0:
                raise PipelineException('found no fasta files in directory "{}"'.format(input_dir))
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
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

        self.complete_step(log, output_dir)
        return output_dir

    def step_06_create_otu_table(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            otus_fp, *_ = glob.glob(os.path.join(input_dir, '*rad3.uchime.fasta'))
            #Get raw reads from earlier steps
            if self.multiple_runs is True:
                input_fps = self.concat_multiple_runs_for_step_06(self.work_dir, output_dir, log)
            elif self.paired_ends is True:
                log.info('Concatenating from step_01_2')
                input_fps = glob.glob(os.path.join(self.work_dir, 'step_01_2*', '*.assembled*.fastq.gz'))
                if self.combine_final_results is True:
                    input_fps = self.concat_all_samples_for_step_06(input_fps, output_dir, log)
            elif self.paired_ends is False and self.cutadapt_min_length != -1:
                log.info('Concatenating from step_01_1')
                input_fps = glob.glob(os.path.join(self.work_dir, 'step_01_1*', '*.fastq.gz'))
                if self.combine_final_results is True:
                    input_fps = self.concat_all_samples_for_step_06(input_fps, output_dir, log)
            elif self.paired_ends is False and self.cutadapt_min_length == -1:
                log.info('Concatenating from step_01_copy_and_compress')
                input_fps = glob.glob(os.path.join(self.work_dir, 'step_01*', '*.fastq.gz'))
                if self.combine_final_results is True:
                    input_fps = self.concat_all_samples_for_step_06(input_fps, output_dir, log)
            if len(input_fps) == 0:
                raise PipelineException('found no .fastq.gz files in directory "{}"'.format(os.path.join(self.work_dir, 'step_02')))
            log.info('input_fps: {}'.format(input_fps))
            for input_fp in input_fps:
                fasta_name = re.sub(
                                string=os.path.basename(input_fp),
                                pattern='\.fastq',
                                repl='.fasta')
                fasta_name = re.sub(
                                string=fasta_name,
                                pattern='\.gz$',
                                repl='')
                fasta_fp = os.path.join(output_dir, fasta_name)
                log.info('convert fastq file\n\t%s\nto fasta file\n\t%s', input_fp, fasta_fp)
                run_cmd([
                    self.vsearch_executable_fp,
                    '--fastq_filter', input_fp,
                    '--fastaout', fasta_fp
                    ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
                otu_table_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq',
                        repl='.uchime.otutab.txt'
                    )
                )
                otu_table_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=otu_table_fp,
                        pattern='\.gz$',
                        repl=''
                    )
                )
                otu_table_biom_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq',
                        repl='.uchime.otutab.biom'
                    )
                )
                otu_table_biom_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=otu_table_biom_fp,
                        pattern='\.gz$',
                        repl=''
                    )
                )
                run_cmd([
                        self.vsearch_executable_fp,
                        '--usearch_global', fasta_fp,
                        '--db', otus_fp,
                        '--id', '0.97',
                        '--biomout', otu_table_biom_fp,
                        '--otutabout', otu_table_fp,
                        '--sizein',
                        '--sizeout',
                        '--threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

                os.remove(fasta_fp)
                if self.multiple_runs is True or self.combine_final_results is True:
                    os.remove(input_fp)

        self.complete_step(log, output_dir)
        return output_dir

    def concat_multiple_runs_for_step_06(self, work_dir, output_dir, log):
        log.info('Concatenating raw reads from multiple runs')
        step_num = ''
        if self.paired_ends is True:
            log.info('Concatenating from step_01_2')
            input_glob = glob.glob(os.path.join(self.work_dir, 'step_01_2*', '*run1*.assembled*.fastq.gz'))
            step_num = '01_2'
        elif self.paired_ends is False and self.cutadapt_min_length != -1:
            log.info('Concatenating from step_01_1')
            input_glob = glob.glob(os.path.join(self.work_dir, 'step_01_1*', '*run1*.fastq.gz'))
            step_num = '01_1'
        elif self.paired_ends is False and self.cutadapt_min_length == -1:
            log.info('Concatenating from step_01_copy_and_compress')
            input_glob = glob.glob(os.path.join(self.work_dir, 'step_01*', '*run1*.fastq.gz'))
            step_num = '01'
        run1_fps = sorted(input_glob)
        log.info('run1 fps: "%s"', str(run1_fps))
        input_fps = []
        for run in run1_fps:
            sample_name = os.path.basename(run).split('_run1')[0]
            sample_glob = os.path.join(work_dir, 'step_%s*' % step_num, '*%s*.fastq.gz*' % sample_name)
            sample_list = sorted(glob.glob(sample_glob))
            log.info('Runs to be concatenated together: "%s"', str(sample_list))
            output_file = os.path.join(output_dir, '%s_concat_runs.fastq.gz' % sample_name)
            with open(output_file, 'wb') as outfile:
                for sample in sample_list:
                    with open(sample, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile) 
            input_fps.append(output_file)
        if self.combine_final_results is True:
            input_fps = self.concat_all_samples_for_step_06(input_fps, output_dir, log)
        return input_fps
        
    def concat_all_samples_for_step_06(self, input_fps, output_dir, log):
        log.info('Concatenating raw reads from all samples')
        log.info('Sample fps for concat_all_samples_for_step_06: "%s"', str(input_fps))
        log.info('unzipping input files')
        uncompressed_input_fps = ungzip_files(*input_fps, target_dir=output_dir)
        if self.multiple_runs is True:
            for i in input_fps:
                os.remove(i)
        log.info('uncompressed file list: "{}"'.format(uncompressed_input_fps))
        combined_fp = '{}.fastq'.format(self.get_combined_file_name(output_dir))
        with open(combined_fp, 'w') as combined_out:
            for input_fp in uncompressed_input_fps:
                input_basename = os.path.basename(input_fp)
                input_name = re.sub(
                            string=input_basename,
                            pattern='_trimmed',
                            repl='')
                input_name = re.sub(
                            string=input_name,
                            pattern='_merged',
                            repl='')
                input_name = re.sub(
                            string=input_name,
                            pattern='\.fastq$',
                            repl='')
                if self.multiple_runs is True:
                    input_name = re.sub(
                                string=input_name,
                                pattern='_concat_runs',
                                repl='')
                with open(input_fp, 'r') as f:
                    for line_ct, l in enumerate(f):
                        if line_ct % 4 == 0:
                            pass
                            l = l[1:]
                            l = '@{}:{}'.format(input_name, l)
                        combined_out.write(l)
                os.remove(input_fp)
        ret_arr = []
        ret_arr.append(combined_fp)
        return ret_arr

    def get_combined_file_name(self, output_dir):
        combined_name = 'total_combined'
        if self.cutadapt_min_length != -1:
            combined_name = '{}_trimmed'.format(combined_name)
        if self.paired_ends is True:
            combined_name = '{}_assembled'.format(combined_name)
        combined_name = '{}_ee{}_trunc{}'.format(combined_name, self.vsearch_filter_maxee, self.vsearch_filter_trunclen)
        combined_fp = os.path.join(output_dir, combined_name)
        return combined_fp
        

if __name__ == '__main__':
    main()
