

# CALCULATE PERFORMANCE METRICS OF CLASSIFICATION MODELS


# Load library to run code in parallel
library(doParallel)

# Load necessary scripts
source('R/packages.R')
source('R/auxiliary_functions.R')
source('R/main_functions.R')

################################################################################################################
# PREPARE DATA

# Select dataset to study (Gandal, Gupta or Wright)
dataset_name = 'Gandal'

# Load data
transcriptomic_data = load_transcriptomic_dataset(dataset_name)
preprocessed_data = transcriptomic_data$preprocessed_data
classification_dataset = transcriptomic_data$classification_dataset
modules_dataset = transcriptomic_data$modules_dataset
biased_classification_model = transcriptomic_data$biased_classification_model
unbiased_classification_model = transcriptomic_data$unbiased_classification_model
SFARI_dataset = read.csv('Results/SFARI_dataset.csv') %>% dplyr::rename('gene-score' = gene.score)

# Add SFARI information to preprocessed_data
genes_info = preprocessed_data$genes_info %>% left_join(SFARI_dataset, by = 'ID') %>%
             mutate(gene.score = ifelse(is.na(`gene-score`) & Neuronal==0, 'Others',
                                        ifelse(is.na(`gene-score`), 'Neuronal', `gene-score`))) %>%
             mutate(Group = factor(ifelse(gene.score %in% c('Neuronal','Others'), gene.score, 'SFARI'),
                                   levels = c('SFARI', 'Neuronal', 'Others'))) %>%
             left_join(preprocessed_data$datGenes %>% mutate('ID' = rownames(.)) %>% 
                       dplyr::select(ID, hgnc_symbol), by = 'ID')
preprocessed_data$genes_info = genes_info

################################################################################################################
# RUN CLASSIFICATION MODEL

# Define number of iterations and seeds
n_iter = 100
seeds = 1:n_iter*1000 # We generate 100 seeds inside the model from the original one, this way they don't overlap 
                      # between iterations

# Setup parallel backend to use many processors
cores = detectCores()
cl = makeCluster(cores-2) # Leaving 2 free cores
registerDoParallel(cl)

# Create matrix with the performance metrics from each iteration
performance_metrics = foreach(i=seeds, .combine=rbind, .errorhandling = 'pass') %dopar% {
  
  source('R/packages.R')
  
  # ###################################################################################
  # # To shuffle SFARI labels uncomment this block
  # set.seed(i)
  # random_classification_dataset = classification_dataset
  # random_classification_dataset$SFARI = sample(classification_dataset$SFARI)
  # random_SFARI_dataset = SFARI_dataset[SFARI_dataset$ID %in% rownames(random_classification_dataset),]
  # random_SFARI_dataset = random_SFARI_dataset %>% mutate(ID = ID %>% as.character)
  # random_SFARI_dataset$ID[!is.na(random_SFARI_dataset$`gene-score`)] =
  #   rownames(random_classification_dataset[random_classification_dataset$SFARI,])
  # random_SFARI_dataset = random_SFARI_dataset %>% distinct(ID, .keep_all = TRUE)
  # ###################################################################################

  # Run model
  results = classification_model(preprocessed_data, classification_dataset, SFARI_dataset, i, FALSE)

  # Save partial results
  temp = as.vector(c('seed' = i, unlist(results$pm_final_model)))
  
}

# Stop clusters
stopCluster(cl)

# Rename columns
colnames(performance_metrics) = c('Seed','Accuracy','Precision','Recall','F1','AUC','MLP','Balanced_Accuracy')

# Set path corresponding to dataset
dataset_path = 'Results/'
if(dataset_name != 'Gandal') dataset_path = paste0('supportingDatasets/', dataset_name, '/' , dataset_path)

# Save results
write.csv(performance_metrics, file = paste0(dataset_path,'performance_biased_model.csv'), row.names=FALSE)
