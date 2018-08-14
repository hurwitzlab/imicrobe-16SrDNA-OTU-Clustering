#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH -p development
#SBATCH -J test-16S
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user mattmiller899@email.arizona.edu

module load tacc-singularity

singularity run $WORK/16S_pipeline/imicrobe-16SrDNA-OTU-Clustering.img \
    --input-dir $WORK/16S_pipeline/data_in/ \
    --pear-min-overlap 10 \
    --pear-max-assembly-length 500 \
    --pear-min-assembly-length 220 \
    --vsearch-filter-maxee 1 \
    --vsearch-filter-trunclen 245 \
    --vsearch-derep-minuniquesize 4 \
    --core-count 4 \
    --multiple-runs \
    --cutadapt-min-length 100

#arg_parser.add_argument('--forward-primer', default='ATTAGAWACCCVNGTAGTCC',
#                            help='forward primer to be clipped by cutadapt')

#    arg_parser.add_argument('--reverse-primer', default='TTACCGCGGCKGCTGGCAC',
#                            help='reverse primer to be clipped by cutadapt')
#${INPUT_DIR} ${OUTPUT_DIR} ${core_count} ${cutadapt_min_length} ${pear_min_overlap} ${pear_max_assembly_length} ${pear_min_assembly_length} ${vsearch_filter_maxee} ${vsearch_filter_trunclen} ${vsearch_derep_minuniquesize} ${forward_primer} ${reverse_primer} ${multiple_runs} ${run_help}
