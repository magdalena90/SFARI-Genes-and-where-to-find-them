
# # WORKFLOW FOR GANDAL'S DATASET

# This script uses drake. Because of some compatibility issues between the packages, it does not perform the 
# enrichment analysis for the top modules. See script workflow.R to run this project without using drake


# Load necessary scripts
source('R/packages.R')
source('R/auxiliary_functions.R')
source('R/main_functions.R')
source('R/plan.R')

# Plot the graph of the workflow
vis_drake_graph(plan, targets_only = TRUE)

# Run project
make(plan)
