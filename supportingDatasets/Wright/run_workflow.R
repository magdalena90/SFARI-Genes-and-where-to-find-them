
# WORKFLOW FOR GUPTA'S DATASET
#setwd('~/PhD/Paper/supportingDatasets/Wright')

# Load necessary scripts
source('../../R/packages.R')
source('../../R/auxiliary_functions.R')
source('../../R/main_functions.R')
source('workflow.R')

# Unload drake (just in case)
detach('package:drake', unload=TRUE)

# Run project
workflow()
