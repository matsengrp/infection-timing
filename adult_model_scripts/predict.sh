#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate infection-timing_py2 
set -eu

PROCESSED_DATA_PATH=$1

cd adult_model_scripts

python predict_time_for_infants.py $PROCESSED_DATA_PATH 

