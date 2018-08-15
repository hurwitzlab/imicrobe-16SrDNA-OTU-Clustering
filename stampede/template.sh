#!/bin/bash

echo "Started $(date)"

echo "INPUT_DIR: ${INPUT_DIR} \nOUTPUT_DIR: ${OUTPUT_DIR} \n UCHIME_REF_DB: ${UCHIME_REF_DB} \ncore_count: ${core_count}\n cutadapt_min_length: ${cutadapt_min_length}\npear_min_overlap: ${pear_min_overlap}\npear_max_assembly_length: ${pear_max_assembly_length}\npear_min_assembly_length: ${pear_min_assembly_length}\nvsearch_filter_maxee: ${vsearch_filter_maxee}\nvsearch_filter_trunclen: ${vsearch_filter_trunclen}\n vsearch_derep_minuniquesize: ${vsearch_derep_minuniquesize}\nforward_primer: ${forward_primer}\nreverse_primer: ${reverse_primer}\n multiple_runs: ${multiple_runs}\n"
sh run.sh ${INPUT_DIR} ${OUTPUT_DIR} ${UCHIME_REF_DB} ${core_count} ${cutadapt_min_length} ${pear_min_overlap} ${pear_max_assembly_length} ${pear_min_assembly_length} ${vsearch_filter_maxee} ${vsearch_filter_trunclen} ${vsearch_derep_minuniquesize} ${forward_primer} ${reverse_primer} ${multiple_runs}

echo "Ended $(date)"
exit 0
