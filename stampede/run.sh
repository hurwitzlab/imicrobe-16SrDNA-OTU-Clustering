#!/bin/bash

module load tacc-singularity

set -u
IMG="/work/05066/imicrobe/singularity/imicrobe-16SrDNA-OTU-Clustering.img"
singularity run "$IMG" "$@" 
echo "Done."
