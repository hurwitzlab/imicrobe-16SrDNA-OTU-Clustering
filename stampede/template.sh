#!/bin/bash

echo "Started $(date)"

sh run.sh ${INPUT_DIR} ${OUTPUT_DIR} ${core_count} ${cutadapt_min_length} ${pear_min_overlap} ${pear_max_assembly_length} ${pear_min_assembly_length} ${vsearch_filter_maxee} ${vsearch_filter_trunclen} ${vsearch_derep_minuniquesize} ${forward_primer} ${reverse_primer} ${multiple_runs} ${run_help}

echo "Ended $(date)"
exit 0
