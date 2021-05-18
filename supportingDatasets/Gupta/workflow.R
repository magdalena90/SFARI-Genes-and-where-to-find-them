
# Complete Project Pipeline

load_and_transform_data = function(){
  
  cat('\n\n\nLoading preprocessed data from PhD_Thesis/Dataset_Gupta/RMarkdowns/2.1.Preprocessing_pipeline...\n')
  
  new_SFARI_dataset = read.csv('../../Results/new_SFARI_dataset.csv') %>% dplyr::rename('gene-score' = gene.score)
  
  load('Results/preprocessed_data.RData')
  genes_info = genes_info %>% left_join(new_SFARI_dataset, by = 'ID') %>%
               mutate(gene.score = ifelse(is.na(`gene-score`) & Neuronal==0, 'Others', 
                                          ifelse(is.na(`gene-score`), 'Neuronal', `gene-score`))) %>%
               mutate(`gene-score` = gene.score) %>% #############################################################
               mutate(Group = factor(ifelse(gene.score %in% c('Neuronal','Others'), gene.score, 'SFARI'), 
                                     levels = c('SFARI', 'Neuronal', 'Others'))) %>%
               left_join(datGenes %>% mutate('ID' = rownames(.)) %>% dplyr::select(ID, hgnc_symbol), by = 'ID')
  datGenes = datGenes %>% mutate('ensembl_gene_id' = rownames(.))
  
  Gupta_dataset = list('datExpr' = datExpr, 'datMeta' = datMeta, 'datGenes' = datGenes, 'genes_info' = genes_info,
                       'dds' = dds)  
}

WGCNA_workflow = function() {
  
  cat('\n\n\nRunning WGCNA...\n')
  
  # Run functions
  modules_dataset = perform_WGCNA(Gupta_dataset)
  top_modules_by_Diagnosis = get_top_modules_by_diagnosis(Gupta_dataset, modules_dataset)
  top_modules_by_SFARI = get_top_modules_by_SFARI(new_SFARI_dataset, modules_dataset)
  
  # Save results
  write.csv(modules_dataset, file = 'Results/modules_dataset.csv', row.names = FALSE)
  #write.csv(top_modules_by_Diagnosis, file = 'Results/top_modules_by_Diagnosis.csv', row.names = FALSE)
  #write.csv(top_modules_by_SFARI, file = 'Results/top_modules_by_SFARI.csv', row.names = FALSE)
  
}

classification_model_workflow = function(){
  
  cat('\n\n\nTraining classification models...\n')
  
  # Run functions
  classification_dataset = classification_model_dataset(Gupta_dataset, modules_dataset, running_on_drake = F)
  biased_classification_model = classification_model(Gupta_dataset, classification_dataset, new_SFARI_dataset, 
                                                     123, correct_bias = F)
  unbiased_classification_model = classification_model(Gupta_dataset, classification_dataset, new_SFARI_dataset, 
                                                       123, correct_bias = T)
  
  # Save Results
  write.csv(classification_dataset, file = 'Results/classification_dataset.csv')
  save(biased_classification_model, file = 'Results/biased_classification_model.RData')
  save(unbiased_classification_model, file = 'Results/unbiased_classification_model.RData')
  
}

enrichment_analysis_workflow = function(){
  
  cat('\n\n\nPerforming Enrichment Analysis of top modules...\n')
  
  top_modules_enrichment = top_modules_EA(top_modules_by_Diagnosis, top_modules_by_SFARI, modules_dataset, classification_dataset)
  save(top_modules_enrichment, file = 'Results/top_modules_enrichment.RData')
  
}

workflow = function(){
  
  load_and_transform_data()
  WGCNA_workflow()
  classification_model_workflow()
  # enrichment_analysis_workflow()
  
  cat('\n\n\nFinished. Results can be found in the Results folder\n')
  
}
