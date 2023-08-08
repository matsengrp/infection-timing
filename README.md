# infection-timing

The goal of this project is to use viral sequence diversity to estimate HIV infection timing for infants using a variety of linear regression methods including Bayesian hierarchical regression and least-absolute-deviation regression. 

# Install
Everything R 4.1.1 based. R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
conda env create -f environment.yml
conda activate infection-timing 
```

If you would like to run analyses involving previously-trained adult-specific least-absolute-deviation models and code (Puller et. al, PLoS Computational Biology 2017), we have built a separate Python-2-based conda environment which can be activated as follows:

```bash 
conda env create -f py2_environment.yml
conda activate infection-timing_py2 
```

# Analysis outline: 

__Table of Contents:__

* [Data preparation](#data_preparation)
* [Model training](#model-training)
    * [Summary of model types](#model-types)
* [Calculating Bayes Factors](#calculating-bayes_factors)
* [Model evaluation](#model-evaluation)
    * [Summary of model evaluation options](#model-evaluation-options)
* [Model validation](#model-validation)
    * [Summary of model validation options](#model-validation-options)
* [Model prediction](#model-prediction)
* [Plot results](#plot-results)
* [Supplementary analyses](#supplementary-analyses)
    * [Model coverage estimates](#model-coverage-estimates)
    * [Model posterior predictive checks](#model-posterior-predictive-checks)

## Data preparation

All of the models presented here leverage HIV viral diversity measures (which we refer to as ADP) to estimate time since infection. The training and testing sequencing data sets used here were originally processed using this [data processing script](preprocess_sequence_data.R). 

### Preparation for making model predictions using your own data with our pre-trained models

We have provided our processed datasets [here](), however, if you would like to make predictions using your own data, make sure your data is in csv format and contains the following required columns:

* `Sample`: string containing patient identifier, gene region identifier, and sequencing run number (i.e. `1_F1_R1` would correspond to individual 1, sequencing fragment/gene region 1, and sequencing run 1)
* `ptnum`: patient identifier
* `pass2_APD`: viral diversity measure
* `Fragment`: gene region identifier (i.e. `1` corresponds to gene region 1, `2` corresponds to gene region 2, etc.)

After formatting your data to include these columns, if you would like to make predictions using our pre-trained models, proceed to the [Model prediction](#model-prediction) section.

### Preparation for training models using your own data

If you would like to train your own models using your own data, you will need several additional columns in your data:

* `month_visit`: the month (post-birth) at which the diversity was sampled
* `vload`: the viral load measured at the specified time point
* `incat_hiv`: whether the individual was known to be infected in-utero or after birth

After formatting your data to include these additional columns, if you would like to train your own models, proceed to the [Model training](#model-training) section.

## Model training

### Model types

Here is a summary of the main model options. You can find a description of additional model type options [here](scripts/stan_models/README.md) and [here](scripts/README.md)

| **Model name (same as manuscript)** | Infant-trained hierarchical model | infant-trained linear models | adult-trained linear models |
|-------------------------------------------------------------------------|
| **Model type** | Bayesian hierarchical regression model | Least absolute deviation linear regression model | Least absolute deviation linear regression model |
| **Data input** | Viral diversity measures (sampled from any of the three gene regions) | Viral diversity measures (sampled from any of the three gene regions), *but a unique model will be trained for each gene region in which data is sampled* | Viral diversity measures (sampled from any of the three gene regions), *but a unique model will be trained for each gene region in which data is sampled* | 
| **Data output** | Distributions of estimated times since infection | Point estimates of time since infection | Point estimates of time since infection |
| **Includes individual- and gene-region-specific slope-modifying terms? | Yes | No | No |
| **Model count** | 1 | 1 per gene region region | 1 per gene region |
| **Script used for model training** | [fit_infant_model](fit_infant_model.sh) | [fit_adult_style_model](fit_adult_style_model.sh) | NA (this model is pre-trained using adult data within Puller et. al, PLoS Comp Bio 2017) |
| **Location where trained models will be stored** | [here](scripts/stan_models/model_fits) | [here](scripts/adult_style_models/model_fits) | NA |

### Workflow

0. Download the training cohort data set using the [link]() provided in our manuscript or use your own data (but make sure the data is formatted the same as ours)
1. Edit the [config](config/config.R) and [file paths](config/file_paths.R) files to be project and/or computer specific. See the [README](config/README.md) for more details.
2. Train model using the model fitting script for your desired model (see table above; for the infant hierarchical model, use [this](fit_infant_model.sh) and for the infant linear models, use [this](fit_adult_style_model.sh). Both of these scripts can be run locally or on a cluster. If using the [fit_infant_model.sh](fit_infant_model.sh) script, it takes a single "time correction type" argument -- options are `uniform`, `beta`, `none`, or `beta_laplace` (see [here](scripts/stan_models/MODEL_DESCRIPTION.md) for more details; the `beta` argument is the one we use throughout the manuscript). The other model training scripts do not have any arguments.

*Model fits will be stored in the location described in the table above (depending on the model type).*

## Calculating Bayes Factors 

If you would like to explore whether APD slopes vary in certain contexts (e.g. by individual, by gene region, etc.), you can use [this script](scripts/bayes_factor_test.R) to calculate Bayes Factors comparing two Bayesian models. This script takes 2 arguments:

1. the time correction type -- options are `uniform`, `beta`, or `none` (e.g. `beta` is the one we use throughout the manuscript) 
2. variation type -- options are described below:

| **Variation type argument** | **Description** |
|-----------------------------------------------|
| `no_fragment` | Compare a model in which APD slopes vary by gene region to a model in which they don't |
| `no_subject` | Compare a model in which APD slopes vary by individual to a model in which they don't |
| `by_infection_time` | Compare a model in which APD slopes vary by mode of infection (i.e. relative timing of infection; either before birth or after birth) to a model in which they don't |
| `with_vload` | Compare a model in which APD slopes vary with viral load to a model in which they don't |
| `with_max_vload` | Compare a model in which APD slopes vary with maximum viral load to a model in which they don't |
| `with_cd4_percent` | Compare a model in which APD slopes vary with percentage of CD4 cells to a model in which they don't |
| `with_cd4_rate` | Compare a model in which APD slopes vary with rate of CD4 cell count decline to a model in which they don't |
| `with_percent_cd4_rate` | Compare a model in which APD slopes vary with rate of CD4 cell percentage decline to a model in which they don't |

This analysis will save a file containing the results in a directory called `bayes_factor_results` located in the indicated `OUTPUT_PATH` as specified in the [config](config) files

## Model evaluation

If you would like to evaluate the performance of a model using various subsets of the training data set (i.e. leave one out cross validation), you can use the [hierarchical model evaluation script](predict_with_infant_model_loocv.sh) or the [linear model evaluation script](predict_with_adult_style_model_loocv.sh). 

Both of these scripts can be run locally or on a cluster. If using the hierarchical model evaluation script, it takes a single "time correction type" argument -- options are `uniform`, `beta`, `none`, or `beta_laplace` (see [here](scripts/stan_models/MODEL_DESCRIPTION.md) for more details; the `beta` argument is the one we use throughout the manuscript). The other model evaluation scripts do not have any arguments.

All output files will be located in a directory called `model_loocv` within the `OUTPUT_PATH` as specified in the [config](config) file

## Model prediction 

If you would like to use the model to make predictions on a new data set and/or validate the model on a testing data set, you can follow these steps:
    
1. Download the processed testing data set, and make sure it contains the same column names as the training data set 
2. Edit [config](config/config.R) file to be project and/or computer specific (be sure to add the path to the file from the last step)
2. Run the [infant-trained hierarchical model prediction script](predict_with_infant_model.sh), [infant-trained linear model prediction script](predict_with_adult_style_model.sh), or the [adult-trained linear model prediction script](predict_with_adult_model.sh). Each of these scripts require different arguments: 

The [infant-trained hierarchical model prediction script](predict_with_infant_model.sh) requires two arguments:

    1. the time correction type--options are `none`, `beta`, or `uniform`
    2. the path to the dataset you want to make predictions for

The [infant-trained linear model prediction script](predict_with_adult_style_model.sh) requires one argument:

    1. the path to the dataset you want to make predictions for

The [adult-trained linear model prediction script](predict_with_adult_style_model.sh) requires two arguments:

    1. the path to the dataset you want to make predictions for
    2. whether the true time since infection is known--options are `TRUE` or `FALSE`

All output files will be located in a directory called `model_predictions` within the `OUTPUT_PATH` as specified in the [config](config) file

## Plot results

Plot [figures](plotting_scripts/manuscript_plots) from the manuscript.

Also, have a look at the plotting [README](plotting_scripts/README.md) for more details.

## Supplementary analyses

### Model coverage estimates

If you would like to compute the coverage interval for the predictions from the infant-trained Bayesian hierarchical regression model, you can run the [coverage script](analysis_scripts/get_coverage.R)

### Model posterior predictive checks

If you would like to conduct a series of posterior predictive checks for the infant-trained Bayesian hierarchical regression model, you can run the [posterior predictive checks script](posterior_pred_check.R)

# About the analysis

With this analysis, we want to quantify rates of viral diversification during infant HIV infection and assess the accuracy of using these diversification rates to estimate time since infection.  

See the manuscript for specific model and methods details: 

Russell, M. L., Fish, C. S., Drescher, S., Cassidy, N. A. J., Chanana, P., Benki-Nugent, S., Slyker, J., Mbori-Ngacha, D., Bosire, R., Richardson, B., Wamalwa, D., Maleche-Obimbo, E., Overbaugh, J., John-Stewart, G., Matsen IV, F. A., Lehman, D. A. (In preparation). Using viral sequence diversity to estimate time of HIV infection in infants.

The following packages were especially helpful in our analyses:

- rstan (Stan Development Team, 2020)
- rstanarm (Goodrich, et. al, 2023)
- data.table (Dowle and Srinivasan, 2021)
- tidyverse (Wickham et. al, 2019) 
- doParallel (Corporation and Steve Weston, 2020)
- cowplot (Wilke, 2020)
