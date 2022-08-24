#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing 
set -eu

DATA_PATH=$1 
TIME_KNOWN=$2

RESULTS_PATH=$(Rscript $PWD/adult_model_scripts/preprocess_infant_data.R $DATA_PATH $TIME_KNOWN)

echo $RESULTS_PATH
