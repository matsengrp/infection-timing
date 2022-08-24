#!/bin/bash

set -eu

DATA_PATH=$1 
TIME_KNOWN=$2

RESULTS_PATH=$(bash adult_model_scripts/preprocess.sh $DATA_PATH $TIME_KNOWN)
bash adult_model_scripts/predict.sh $RESULTS_PATH

echo "done"
