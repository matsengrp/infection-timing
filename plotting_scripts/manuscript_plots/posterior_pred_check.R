library(rstan)
library(rstanarm)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)

TIME_CORRECTION_TYPE <<- 'beta'
NCPU <<- 2
PREPROCESS_DATA <<-  TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)
data$fragment_id = data$fragment

model = load_model_fit()
post = rstan::extract(model)
transformed_posteriors = transform_posterior_matrix_to_dataframe(data, post$observed_time_since_infection_rep, observed_time = TRUE)
raw = unique(transformed_posteriors[, -c('iteration', 'time_since_infection_posterior_draw')])

for (stat in c('median', 'mean', 'sd', 'max', 'min')){
    print(calculate_posterior_pred_check(get(stat), transformed_posteriors, raw))
    assign(paste0(stat, '_plot'), plot_posterior_pred_check(get(stat), stat, transformed_posteriors, raw, write_plot = FALSE))
}

all = align_plots(plot_median, plot_sd, plot_max, plot_min, align = 'vh', axis = 'lr')

grid = plot_grid(all[[1]], all[[2]], all[[3]], all[[4]], nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 35) 

name = paste0('plots/manuscript_figs/all_post_pred_checks.pdf')
ggsave(name, grid, width = 25, height = 20, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


