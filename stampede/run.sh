#!/bin/bash

module load tacc-singularity
module load launcher

set -u

function ADVANCED_USAGE() {
    echo "Usage: pipeline.py [-h] -i INPUT_DIR -w WORK_DIR"
    echo "      -b BARCODE_LENGTH -m MAPPING_FILE [-p PAIRED_DIR]"
    echo
    echo "Required arguments:"
    echo "  -i INPUT_DIR            path to the input directory (or single file)"
    echo "  -w WORK_DIR         path to the output directory"
    echo "  -b BARCODE_LENGTH       length of the barcodes (int)"
    echo "  -m MAPPING_FILE         path to file containing information about the barcodes. Must be in the QIIME mapping file format"
    echo
    echo "Optional arguments:"
    echo "  -h              show this help message and exit"
    echo "  -p PAIRED_DIR              path to the paired end directory (or single file)"
    echo
    exit 1
}


function USAGE() {
    echo "Usage: pipeline.py [-h] -i INPUT_DIR -w WORK_DIR"
    echo "      -b BARCODE_LENGTH -m MAPPING_FILE [-p PAIRED_DIR]"
    echo "Required arguments:"
    echo "  -b BARCODE_LENGTH"
    echo "  -m MAPPING_FILE"
    echo "  -i INPUT_DIR"
    echo "  -w WORK_DIR"
    echo
    echo "Options:"
    echo "  -h"
    echo "  -p PAIRED_DIR"
    echo
    exit 1
}
INPUT_DIR=$1
OUTPUT_DIR=$2
core_count=$3
cutadapt_min_length=$4
pear_min_overlap=$5
pear_max_assembly_length=$6
pear_min_assembly_length=$7
vsearch_filter_maxee=$8
vsearch_filter_trunclen=$9
vsearch_derep_minuniquesize=$10
forward_primer=$11
reverse_primer=$12
multiple_runs=$13
run_help=$14

if [ $run_help = true ]; then
    ADVANCED_USAGE
fi

if [[ $INPUT_DIR -eq "this is needed" ]]; then
    echo "INPUT_DIR is required"
    exit 1
fi

if [[ $OUTPUT_DIR -eq "this is needed" ]]; then
    echo "OUTPUT_DIR is required"
    exit 1
fi




echo "starting directory : `pwd`"
echo "`ls -l`"
echo "input directory    : ${INPUT_DIR}"
echo "output directory   : ${OUTPUT_DIR}"

export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_WORKDIR=${OUTPUT_DIR}
export LAUNCHER_RMI=SLURM

export LAUNCHER_JOB_FILE=`pwd`/launcher_jobfile_${SLURM_JOB_ID}
echo ${LAUNCHER_JOB_FILE}

xz --decompress imicrobe-16SrDNA-OTU-Clustering.img.xz

echo "`ls -l`"

if [ $multiple_runs = true ]; then
    time singularity run imicrobe-16SrDNA-OTU-Clustering.img \
      --input-dir ${INPUT_DIR} \
      --work-dir ${OUTPUT_DIR} \
      --core-count ${core_count} \
      --cutadapt-min-length ${cutadapt_min_length} \
      --pear-min-overlap ${pear_min_overlap} \
      --pear-max-assembly-length ${pear_max_assembly_length} \
      --pear-min-assembly-length ${pear_min_assembly_length} \
      --vsearch-filter-maxee ${vsearch_filter_maxee} \
      --vsearch-filter-trunclen ${vsearch_filter_trunclen} \
      --vsearch-derep-minuniquesize ${vsearch_derep_minuniquesize} \
      --forward-primer ${forward_primer} \
      --reverse-primer ${reverse_primer} \
      --multiple-runs
else
    time singularity run imicrobe-16SrDNA-OTU-Clustering.img \
      --input-dir ${INPUT_DIR} \
      --work-dir ${OUTPUT_DIR} \
      --core-count ${core_count} \
      --cutadapt-min-length ${cutadapt_min_length} \
      --pear-min-overlap ${pear_min_overlap} \
      --pear-max-assembly-length ${pear_max_assembly_length} \
      --pear-min-assembly-length ${pear_min_assembly_length} \
      --vsearch-filter-maxee ${vsearch_filter_maxee} \
      --vsearch-filter-trunclen ${vsearch_filter_trunclen} \
      --vsearch-derep-minuniquesize ${vsearch_derep_minuniquesize} \
      --forward-primer ${forward_primer} \
      --reverse-primer ${reverse_primer} 
fi

#sleep 10
#export LAUNCHER_PPN=2

#$LAUNCHER_DIR/paramrun
echo "Ended launcher"
