# Configuration Files

These files contain variables for repository path, output data path, data paths, etc. 
The variables must be changed to be computer/project specific. 

Specifically, the following change is required within the [config.R](config.R) file:

* `PROJECT_PATH` (location of the `infection-timing` repository) and `OUTPUT_PATH` (location where output files should be stored) variables

Within the [file_paths.R](file_paths.R) file, the following changes are required:

* `TRAINING_INFANT_DATA_PATH` (file path for the training data set) and `TRAINING_INFANT_METADATA_PATH` (file path for the training data set meta data)
* `TESTING_INFANT_DATA_PATH` (file path for the testing data set) and `TESTING_INFANT_DATA_TRUE_TIME_PATH` (file path for the true infection times for the training data set)
