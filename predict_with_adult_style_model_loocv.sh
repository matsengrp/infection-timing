#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing
set -eu

Rscript scripts/cross_validate_adult_style_model.R 

echo "done"
