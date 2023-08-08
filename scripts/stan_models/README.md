Currently, there are several stan models: 

* a model [without a time correction term](multilevel_infant_model.stan)
* a model with a [uniform time correction term](multilevel_infant_model_time_error_uniform.stan)
* a model with a [beta distribution time correction term](multilevel_infant_model_time_error_beta.stan)

These models are identical otherwise. The model that we refer to as the infant-trained hierarchical model is the model with a beta distribution time correction term.
For a description of this model and its parameters, see the model description [here](MODEL_DESCRIPTION.md)

The errors for the above model with a beta distribution time correction term are normally distributed. We also have a similar model in which errors are Laplace distributed [here](multilevel_infant_model_time_error_beta_laplace.stan).

We also have several models that include (or lack) various additional slope-varying effects. In our analyses, we compare the model with a beta distribution time correction term and normally distributed errors ([here](multilevel_infant_model_time_error_beta.stan)) to alternative models by calculating Bayes factors. This original model contains individual-specific and gene-region-specific slope modifying terms. We can compare this model to the following alternative models:

* a model [lacking an individual-specific slope modifying term](multilevel_infant_model_time_error_beta_no_subject_variation.stan)
* a model [lacking a gene-region-speicific slope modifying term](multilevel_infant_model_time_error_beta_no_fragment_variation.stan)
* a model [containing an additional mode-of-infection-specific slope modifying term](multilevel_infant_model_time_error_beta_with_infection_time_variation.stan)
* a model [containing an additional viral load slope modifier](multilevel_infant_model_time_error_beta_with_vload_variation.stan)
* a model [containing an additional CD4 decline rate slope modifier](multilevel_infant_model_time_error_beta_with_percent_cd4_rate_variation.stan)

In addition to these models, we also have two prediction models that allow us to make predictions on new data using a "frozen" model:

* [standard prediction model](predict.stan)
* [single fragment prediction model](predict_one_frag.stan)
