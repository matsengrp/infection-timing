#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing 
set -eu

TIME_CORRECTION_TYPE=$1

Rscript scripts/fit_model.R $TIME_CORRECTION_TYPE

echo "done"
