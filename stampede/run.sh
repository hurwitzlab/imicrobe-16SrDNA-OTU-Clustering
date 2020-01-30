#!/bin/bash

set -u

IMG="/work/05066/imicrobe/singularity/imicrobe-16SrDNA-OTU-Clustering-0.1.0.img"
singularity run "$IMG" "$@" 
echo "Done."
