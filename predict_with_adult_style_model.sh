#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing
set -eu

DATA_PATH=$1

Rscript scripts/predict_adult_style_model.R $DATA_PATH

echo "done"
