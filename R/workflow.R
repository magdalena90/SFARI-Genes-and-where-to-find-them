
# Complete Project Pipeline

preprocessing_workflow = function(){
  
  cat('\n\n\nPreprocessing input data...\n')
  
  # Run functions
  SFARI_dataset = preprocess_SFARI(version='new')
  NCBI_dataset = preprocess_NCBI()
  GO_neuronal_dataset = preprocess_GO_annotations()
  preprocessed_dataset = preprocess_Gandal(GO_neuronal_dataset, NCBI_dataset, SFARI_dataset)
  
  # Save results
  write.csv(SFARI_dataset, file = 'Results/SFARI_dataset.csv', row.names = FALSE)
  write.csv(NCBI_dataset, file = 'Results/NCBI_dataset.csv', row.names = FALSE)
  write.csv(GO_neuronal_dataset, file = 'Results/GO_neuronal_dataset.csv', row.names = FALSE)
  save(preprocessed_dataset, file = 'Results/preprocessed_data.RData')
  
}

WGCNA_workflow = function() {
  
  cat('\n\n\nRunning WGCNA...\n')
  
  # Load input
  load('Results/Gandal_dataset.RData')
  SFARI_dataset = read.csv('Results/SFARI_dataset.csv') %>% dplyr::rename('gene-score' = gene.score)
  
  # Run functions
  modules_dataset = perform_WGCNA(Gandal_dataset)
  top_modules_by_Diagnosis = get_top_modules_by_diagnosis(Gandal_dataset, modules_dataset)
  top_modules_by_SFARI = get_top_modules_by_SFARI(SFARI_dataset, modules_dataset)
  
  # Save results
  write.csv(modules_dataset, file = 'Results/modules_dataset.csv', row.names = FALSE)
  write.csv(top_modules_by_Diagnosis, file = 'Results/top_modules_by_Diagnosis.csv', row.names = FALSE)
  write.csv(top_modules_by_SFARI, file = 'Results/top_modules_by_SFARI.csv', row.names = FALSE)
  
}

classification_model_workflow = function(){
  
  cat('\n\n\nTraining classification models...\n')
  
  # Load input
  load('Results/Gandal_dataset.RData')
  SFARI_dataset = read.csv('Results/SFARI_dataset.csv') %>% dplyr::rename('gene-score' = gene.score)
  modules_dataset = read.csv('Results/modules_dataset.csv')
  
  # Run functions
  classification_dataset = classification_model_dataset(Gandal_dataset, modules_dataset, running_on_drake = F)
  biased_classification_model = classification_model(Gandal_dataset, classification_dataset, SFARI_dataset, 
                                                     123, correct_bias = F)
  unbiased_classification_model = classification_model(Gandal_dataset, classification_dataset, SFARI_dataset, 
                                                       123, correct_bias = T)
  
  # Save Results
  write.csv(classification_dataset, file = 'Results/classification_dataset.csv')
  save(biased_classification_model, file = 'Results/biased_classification_model.RData')
  save(unbiased_classification_model, file = 'Results/unbiased_classification_model.RData')
  
}

enrichment_analysis_workflow = function(){
  
  cat('\n\n\nPerforming Enrichment Analysis of top modules...\n')
  
  top_modules_enrichment = top_modules_EA(top_modules_by_Diagnosis, top_modules_by_SFARI, modules_dataset, 
                                          classification_dataset)
  save(top_modules_enrichment, file = 'Results/top_modules_enrichment.RData')
  
}

workflow = function(){

  preprocessing_workflow()
  WGCNA_workflow()
  classification_model_workflow()
  #enrichment_analysis_workflow()
  
  cat('\n\n\nFinished. Results can be found in the Results folder\n')
  
}
