#!/bin/bash

module load tacc-singularity

set -u
IMG="/work/05286/mattmill/stampede2/16S_pipeline/imicrobe-16SrDNA-OTU-Clustering-0.0.2/stampede/imicrobe-16SrDNA-OTU-Clustering.img"
singularity run "$IMG" "$@" 
#sleep 10
#export LAUNCHER_PPN=2

#$LAUNCHER_DIR/paramrun
echo "Ended launcher"
