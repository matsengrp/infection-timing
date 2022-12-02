library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/scripts/bayes_factor_functions.R'))

bf_filepath = get_bayes_factor_results_path()
bf = fread(bf_filepath)
