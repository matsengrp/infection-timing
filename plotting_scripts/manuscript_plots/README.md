# Manuscript figures

This directory contains a script to produce each manuscript figure.
The scripts map to manuscript figures as follows:

- [Figure 1A, Figure 2, Figure S1](plot_apd_time_data_credible_intervals.R): figure comparing APD to observed time (Figure 1A) and overlaid with inferred model APD slopes (Figure 2); also a figure showing credible intervals for inferred APD slopes (Figure S1)
- [Figure 3](plot_vload_rate.R): figure comparing median APD slopes with set-point viral load
- [Figure 5](training_side_by_side.R): figure showing mean absolute error by sequencing region, comparison of true time since infection and model-derived time since infection, and mean residual versus true time since infection
- [Figure 7](testing_side_by_side.R): figure showing mean absolute error by sequencing region using the testing dataset
- [Figure S2](plot_cd4_rate.R): figure comparing median APD slopes with rate of CD4 decline
- [Figure S3](posterior_pred_check.R): figure showing histograms for various model posterior predictive checks
- [Figure S4](training_compare_laplace.R): comparison of the model prediction errors for the training data set when drawing observations from a Laplace distribution and a Normal distribution within the infant-trained hierarchical model formulation 
- [Figure S5 and Figure S6](compare_training_testing_data.R): distribution of APD measures for the three gene-regions for each data set
