
# Drake Plan

plan <- drake_plan (
  
  # Preprocess input data
  SFARI_dataset = preprocess_SFARI(version='new'),
  NCBI_dataset = preprocess_NCBI(),
  GO_neuronal_dataset = preprocess_GO_annotations(),
  Gandal_dataset = preprocess_Gandal(GO_neuronal_dataset, NCBI_dataset, SFARI_dataset),
  
  # WGCNA
  modules_dataset = perform_WGCNA(Gandal_dataset),
  top_modules_by_Diagnosis = get_top_modules_by_diagnosis(Gandal_dataset, modules_dataset),
  top_modules_by_SFARI = get_top_modules_by_SFARI(SFARI_dataset, modules_dataset),
  #top_modules_enrichment = top_modules_EA(top_modules_by_Diagnosis, top_modules_by_SFARI, modules_dataset, 
  #                         classification_dataset),
  
  # Classification Model
  classification_dataset = classification_model_dataset(Gandal_dataset, modules_dataset, running_on_drake = T),
  biased_classification_model = classification_model(Gandal_dataset, classification_dataset, SFARI_dataset, 123, 
                                                     correct_bias = F),
  unbiased_classification_model = classification_model(Gandal_dataset, classification_dataset, SFARI_dataset, 123, 
                                                       correct_bias = T)
  
)

