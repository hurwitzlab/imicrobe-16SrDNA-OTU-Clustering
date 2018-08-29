#!/bin/bash

echo "Started $(date)"

echo "INPUT_DIR: ${INPUT_DIR}"
echo "UCHIME_REF_DB: ${UCHIME_REF_DB}"
echo "core_count: ${core_count}"
echo "cutadapt_min_length: ${cutadapt_min_length}"
echo "pear_min_overlap: ${pear_min_overlap}"
echo "pear_max_assembly_length: ${pear_max_assembly_length}"
echo "pear_min_assembly_length: ${pear_min_assembly_length}"
echo "vsearch_filter_maxee: ${vsearch_filter_maxee}"
echo "vsearch_filter_trunclen: ${vsearch_filter_trunclen}"
echo "vsearch_derep_minuniquesize: ${vsearch_derep_minuniquesize}"
echo "forward_primer: ${forward_primer}"
echo "reverse_primer: ${reverse_primer}"
echo "multiple_runs: ${multiple_runs}"
echo "paired_ends: ${paired_ends}"
echo "debug: ${debug}"
echo "steps: ${steps}"

sh run.sh ${INPUT_DIR} ${UCHIME_REF_DB} ${core_count} ${cutadapt_min_length} ${pear_min_overlap} ${pear_max_assembly_length} ${pear_min_assembly_length} ${vsearch_filter_maxee} ${vsearch_filter_trunclen} ${vsearch_derep_minuniquesize} ${forward_primer} ${reverse_primer} ${multiple_runs} ${paired_ends} ${steps} ${debug}

echo "Ended $(date)"
