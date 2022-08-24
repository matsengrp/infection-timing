#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing 
set -eu

TIME_CORRECTION_TYPE=$1

Rscript scripts/cross_validate.R $TIME_CORRECTION_TYPE

echo "done"
