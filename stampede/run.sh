#!/bin/bash

module load tacc-singularity

set -u
IMG="/work/05286/mattmill/stampede2/16S_pipeline/imicrobe-16SrDNA-OTU-Clustering.img"
singularity run "$IMG" "$@" 
echo "Done."
