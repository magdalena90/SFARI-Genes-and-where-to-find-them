
# WORKFLOW FOR GANDAL'S DATASET

# This script was created as an alternative to the make.R file, which runs using drake


# Load necessary scripts
source('R/packages.R')
source('R/auxiliary_functions.R')
source('R/main_functions.R')
source('R/workflow.R')

# Unload drake
detach('package:drake', unload=TRUE)

# Run project
workflow()
