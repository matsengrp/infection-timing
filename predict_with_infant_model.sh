#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing
set -eu

TIME_CORRECTION_TYPE=$1
DATA_PATH=$2

Rscript scripts/predict.R $TIME_CORRECTION_TYPE $DATA_PATH

echo "done"
