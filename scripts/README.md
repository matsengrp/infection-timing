# Analysis scripts

This directory contains files for running the all of the analyses: 

The main scripts are contained in the following files:

## General functions

* [Model training functions](model_fitting_functions.R) which are sourced by nearly all top-level modeling scripts
* [Model cross validation functions](model_cross_validation_functions.R) which are sourced by nearly all top-level analysis scripts doing cross validation
* [Model prediction functions](model_prediction_evaluation_functions.R) which are sourced by nearly all top-level analysis scripts doing model prediction or validation
* [Linear regression style model fuctions](adult_style_model_functions.R) which are sourced by nearly all top-level modeling scripts involving the linear regression style model

## Model training scripts

* [Model training script](fit_model.R) which is sourced by top level model training script (e.g. ../fit_infant_model.sh)
* [Model training script (for MAP fitting)](fit_model_map.R) 
* [Model training script (for linear regression style model)](fit_adult_style_model.R) which is sourced by top level model training script (e.g. ../fit_adult_style_model.sh)

## Model cross validation scripts

* [Model cross validation script](cross_validate.R) which is sourced by top level model cross validation script (e.g. ../predict_with_infant_model_loocv.sh)
* [Model cross validation script (for MAP fitting)](cross_validate_map.R) which is sourced by top level model cross validation script (e.g. ../predict_with_infant_model_map_loocv.sh)
* [Model cross validation script (for linear regression style model)](cross_validate_adult_style_model.R) which is sourced by top level model cross validation script (e.g. ../predict_with_adult_style_model_loocv.sh)

## Model prediction scripts

* [Model prediction script](predict.R) which is sourced by top level model prediction script (e.g. ../predict_with_infant_model.sh)
* [Model prediction script (for MAP trained model)](predict_map.R)
* [Model prediction script (for linear regression style model)](predict_adult_style_model.R) which is sourced by top level model prediction script (e.g. ../predict_with_adult_style_model.sh)

## Bayes Factor test scripts and functions

* [Bayes factor test functions](bayes_factor_functions.R) which is sourced by the Bayes factor test script
* [Bayes factor test script](bayes_factor_test.R)


Additional analysis scripts are contained in the [analysis scripts](analysis_scripts) directory. See [README](analysis_scripts/README.md) for more details.

Argument specific functions are located in the other directories housed here:

* [Stan models](stan_models) 
* [Model time correction type](model_time_correction_type)

"Frozen" model fits:

* [Infant-trained Bayesian hierarchical models](stan_models/model_fits)
* [Infant-trained linear regression models](adult_style_models/model_fits)
