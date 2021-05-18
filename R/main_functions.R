
# MAIN FUNCTIONS

preprocess_NCBI = function(){
  
  # Data downladed from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
  
  # Load data and filtering tax_id=9606 to keep the genes belonging to the human genome only
  gene2ensembl = fread(file='InputData/NCBI_gene2ensembl_20_02_07.gz') %>% filter(`#tax_id` == 9606)
  gene_info = fread(file='InputData/NCBI_gene_info_20_02_07.gz') %>% filter(`#tax_id` == 9606)
  
  # Merge datasets, rename columns and modify gene_biotype values
  NCBI_dataset = gene2ensembl %>% 
                 dplyr::select(c(2,which(!colnames(gene2ensembl) %in% colnames(gene_info)))) %>% 
                 left_join(gene_info, by='GeneID') %>% 
                 dplyr::select(Ensembl_gene_identifier, Symbol, type_of_gene) %>%
                 distinct(Ensembl_gene_identifier, Symbol, type_of_gene) %>%
                 dplyr::rename('ensembl_gene_id'=Ensembl_gene_identifier, 'gene_biotype'=type_of_gene, 
                               'hgnc_symbol'=Symbol) %>% 
                 mutate(gene_biotype = ifelse(gene_biotype=='protein-coding','protein_coding', gene_biotype))
  
  return(NCBI_dataset)
  
}

preprocess_SFARI = function(version='new'){
  
  ##############################################################################################################
  # LOAD DATASET: Downloaded from https://gene.sfari.org/database/gene-scoring/
  
  if(version=='new'){
    SFARI_dataset = read_csv('InputData/SFARI_genes_01-03-2020.csv')
  } else {
    SFARI_dataset = read_csv('InputData/SFARI_genes_08-29-2019.csv')
  }
  
  ##############################################################################################################
  # ADD ENSEMBL IDs: Complementing the current archive results with the feb2014 archive
  
  # Current version (the latest version when preparing this code was jan2020):
  getinfo = c('hgnc_symbol','ensembl_gene_id','gene_biotype')
  mart = umart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl',
                         host = 'jan2020.archive.ensembl.org')
  datGenes_now = getBM(attributes=getinfo, filters='hgnc_symbol', values=SFARI_dataset$`gene-symbol`, mart=mart)
  
  # feb2014  archive
  mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', 
                 host='feb2014.archive.ensembl.org')
  datGenes_2014 = getBM(attributes=getinfo, filters='hgnc_symbol',values=SFARI_dataset$`gene-symbol`,mart=mart)
  datGenes_2014 = datGenes_2014 %>% filter(!ensembl_gene_id %in% datGenes_now$ensembl_gene_id)
  
  # Merge results
  datGenes = rbind(datGenes_now, datGenes_2014) %>% data.frame %>% arrange(hgnc_symbol)
  
  # Merge SFARI_dataset with Ensembl IDs
  SFARI_dataset = SFARI_dataset %>% 
                  left_join(datGenes, by=c('gene-symbol'='hgnc_symbol')) %>% 
                  mutate(ID = ensembl_gene_id)
  
  
  ##############################################################################################################
  # IMPROVING BIOMART RESULTS
  
  # Filter IDs with 'LRG' format
  SFARI_dataset = SFARI_dataset %>% filter(!grepl('LRG_',ensembl_gene_id))
  
  # Adding biotype manually for genes where BioMart return no results
  SFARI_dataset = SFARI_dataset %>% 
                  mutate(ID = case_when(`gene-symbol`=='MSNP1AS'        ~ 'ENSG00000251593',
                                        `gene-symbol`=='RP11-1407O15.2' ~ 'ENSG00000174093',
                                         TRUE ~ ID),
                         gene_biotype = case_when(`gene-symbol`=='MSNP1AS'        ~ 'processed_pseudogene', 
                                                  `gene-symbol`=='RP11-1407O15.2' ~ 'protein_coding',
                                                   TRUE ~ gene_biotype)) %>%
                  dplyr::select(-ensembl_gene_id)
  
  return(SFARI_dataset)

}

preprocess_GO_annotations = function(){
  
  # Downloaded from http://geneontology.org/
  GO_annotations = read.csv('InputData/genes_GO_annotations.csv')
  
  GO_neuronal_dataset = GO_annotations %>% 
                        filter(grepl('neuro', go_term)) %>% 
                        mutate('ID'=as.character(ensembl_gene_id)) %>% 
                        dplyr::select(-ensembl_gene_id) %>% distinct(ID) %>%
                        mutate('Neuronal' = 1)
  
  return(GO_neuronal_dataset)
  
}

preprocess_Gandal = function(GO_neuronal_dataset, NCBI_dataset, SFARI_dataset){
  
  # 1. Load and prepare data
  datExpr = read.csv('InputData/RNAseq_ASD_datExpr.csv', row.names=1)
  datMeta = read.csv('InputData/RNAseq_ASD_datMeta.csv') %>% 
            mutate(Brain_Region = as.factor(Region)) %>% 
            mutate(Brain_lobe = case_when(Brain_Region %in% c('BA4_6', 'BA9', 'BA24', 'BA44_45') ~ 'Frontal',
                                          Brain_Region %in% c('BA3_1_2_5', 'BA7') ~ 'Parietal',
                                          Brain_Region %in% c('BA38','BA39_40','BA20_37','BA41_42_22') ~ 'Temporal',
                                          TRUE ~ 'Occipital'),
                   Batch = as.factor(gsub('/', '.', RNAExtractionBatch)),
                   Diagnosis = factor(Diagnosis_, levels=c('CTL','ASD'))) %>% 
            dplyr::select(-c(Diagnosis_, Subject_ID_2, RNAExtractionBatch, contains('Picard')))
  
  # 2. Anotate genes with BioMart information
  datGenes = annotate_genes(datExpr, NCBI_dataset)
  
  # 3. Filtering Genes and Samples
  filtered_data = filter_genes_and_samples(datExpr, datGenes, datMeta, threshold = 75)
  
  # 4. Look for unknown sources of batch effects with SVA
  datMeta = perform_sva(filtered_data)
  datExpr = filtered_data[['datExpr']]
  datGenes = filtered_data[['datGenes']]
  rm(filtered_data)
  
  # 5. Normalisation and Differential Expression Analysis
  norm_data = normalisation_and_DEA(datExpr, datGenes, datMeta)
  
  # 6. Correct Batch Effects
  datExpr = correct_batch_effects(norm_data)
  
  # 7. Extract final versions of our datasets
  datExpr = datExpr %>% data.frame
  datMeta = norm_data[['datMeta']]
  DE_info = norm_data[['DE_info']] %>% data.frame
  dds = norm_data$dds
  rm(norm_data)
  
  # 8. Incorporate DEA, Neuronal and SFARI Genes info into a single data frame
  genes_info = DE_info %>% 
               data.frame %>%
               mutate(ID = rownames(.)) %>% 
               left_join(SFARI_dataset, by = 'ID') %>% 
               mutate(`gene-score` = ifelse(is.na(`gene-score`), 'Others', `gene-score`)) %>%
               distinct(ID, .keep_all = TRUE) %>% 
               left_join(GO_neuronal_dataset, by = 'ID') %>% 
               mutate(Neuronal = ifelse(is.na(Neuronal), 0, Neuronal)) %>% 
               mutate(gene.score = ifelse(`gene-score`=='Others' & Neuronal==1, 'Neuronal', `gene-score`), 
                      significant = padj<0.05 & !is.na(padj)) %>%
               left_join(datGenes %>% dplyr::select(-gene_biotype), by = c('ID'='ensembl_gene_id')) %>%
               dplyr::select(-c(`ensembl-id`,`genetic-category`,`gene-symbol`)) %>%
               dplyr::select(c(ID, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, 
                               significant, syndromic, Neuronal ,`gene-score`, gene.score))
  
  
  norm_data = list('datExpr'=datExpr, 'datGenes'=datGenes, 'datMeta'=datMeta, 'genes_info'=genes_info, 'dds'=dds)
  
  return(norm_data)
  
}

perform_WGCNA = function(dataset){
  
  get_mod_colours = function(mods){
    
    n = length(unique(mods))-1
    set.seed(123) ; rand_order = sample(1:n)
    mod_colors = c('gray',gg_colour_hue(n)[rand_order])
    names(mod_colors) = mods %>% table %>% names
    
    return(mod_colors)
  }
  
  # Extract information from input
  datExpr = dataset[['datExpr']]
  datGenes = dataset[['datGenes']]
  datMeta = dataset[['datMeta']]
  
  # Scale-Free Topology
  allowWGCNAThreads()
  best_power = datExpr %>% t %>% pickSoftThreshold(powerVector = seq(1,6,0.1), RsquaredCut=0.8)
  S_sft = datExpr %>% t %>% adjacency(type='signed hybrid', power=best_power$powerEstimate, corFnc='bicor')
  
  # Build dissimilarity matrix from TOM
  dissTOM = S_sft %>% TOMdist
  rownames(dissTOM) = rownames(S_sft)
  colnames(dissTOM) = colnames(S_sft)
  
  # Clustering
  dend = dissTOM %>% as.dist %>% hclust(method='average')
  modules = cutreeDynamic(dend, minClusterSize = 10, distM = dissTOM)
  module_colors = get_mod_colours(modules)
  modules_dataset = data.frame('ID' = rownames(datExpr), 'Module' = module_colors[as.character(modules)])
  
  return(modules_dataset)
  
}

get_top_modules_by_diagnosis = function(dataset, modules_dataset){
  
  # Extract information from input
  datExpr = dataset[['datExpr']]
  datMeta = dataset[['datMeta']]
  
  # Calculate association to diagnosis and extract modules with a correlation magnitude higher than 0.9
  top_modules = module_diagonsis_association(datExpr, datMeta, modules_dataset) %>% 
                arrange(desc(MTcor)) %>% 
                filter(abs(MTcor)>0.9) %>% 
                dplyr::select(Module, MTcor)
  
  return(top_modules)
  
}

get_top_modules_by_SFARI = function(SFARI_dataset, modules_dataset){

  # Prepare input:
  
  EA_dataset = prepare_dataset_for_EA(modules_dataset)
  
  term2gene = EA_dataset %>% 
              left_join(SFARI_dataset %>% dplyr::select(ID, `gene-score`), by = 'ID') %>% 
              mutate('SFARI' = ifelse(is.na(`gene-score`), 'Others', 'SFARI'),
                     `gene-score` = ifelse(!is.na(`gene-score`), paste0('Score ',`gene-score`), 'Others')) %>%
              dplyr::select(-ID, -Module) %>%
              melt(id.vars = 'entrezgene') %>%
              dplyr::rename('term' = value, 'gene' = entrezgene) %>% 
              dplyr::select(term, gene) %>%
              distinct
  
  modules = EA_dataset %>% pull(Module) %>% unique %>% as.character
  
  universe = EA_dataset %>% pull(entrezgene) %>% as.character

  
  # Perform Enrichment:
  
  top_modules = c()
  for(module in modules){
    
    genes_in_module = EA_dataset %>% filter(Module == module) %>% pull(entrezgene) %>% as.character
    
    ORA = enricher(gene = genes_in_module, universe = universe, pAdjustMethod = 'bonferroni', 
                   TERM2GENE = term2gene, pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize = 1000)
    
    if(!is.null(ORA) && ORA@result$p.adjust[ORA@result$ID=='SFARI']<0.01) {
      new_row = c(module, 1-ORA@result$pvalue[ORA@result$ID=='SFARI'], ORA@result$p.adjust[ORA@result$ID=='SFARI'])
      top_modules = top_modules %>% rbind(new_row)
    }
  }
  
  top_modules = top_modules %>% 
                data.frame(stringsAsFactors = FALSE) %>%
                dplyr::rename('Module' = X1, 'Enrichment' = X2, 'adj.pval' = X3) %>%
                mutate(Enrichment = Enrichment %>% as.numeric, adj.pval = adj.pval %>% as.numeric) %>%
                arrange(desc(Enrichment))
                
  
  return(top_modules)

}
  
top_modules_EA = function(top_modules_by_Diagnosis, top_modules_by_SFARI, modules_dataset, classification_dataset){
  
  top_modules_enrichment = list()
  
  # Get top modules
  top_modules = c(top_modules_by_Diagnosis %>% pull(Module) %>% as.character, 
                  top_modules_by_SFARI %>% pull(Module) %>% as.character) %>% 
                unique
  
  # Perform Enrichment Analysis
  EA_dataset = prepare_dataset_for_EA(modules_dataset)
  GSEA_enrichment = enrichment_with_GSEA(EA_dataset, top_modules, classification_dataset, nPerm = 1e5)
  ORA_enrichment = enrichment_with_ORA(EA_dataset, top_modules)
  
  # Get sahred enrichment for each module
  for(module in top_modules){
    
    module_enrichment = list()
    GSEA_enrichment_for_module = GSEA_enrichment[[module]]
    ORA_enrichment_for_module = ORA_enrichment[[module]]
    
    for(dataset in c('KEGG', 'Reactome', 'GO', 'DO', 'DGN')){
      
      GSEA_enrichment_dataset = GSEA_enrichment_for_module[[dataset]] %>% 
                                data.frame %>%
                                dplyr::rename('pvalue_GSEA' = pvalue, 
                                              'p.adjust_GSEA' = p.adjust, 
                                              'qvalues_GSEA' = qvalues)
                              
      ORA_enrichment_dataset = ORA_enrichment_for_module[[dataset]] %>% 
                               data.frame %>%
                               dplyr::rename('pvalue_ORA' = pvalue, 
                                              'p.adjust_ORA' = p.adjust, 
                                              'qvalue_ORA' = qvalue)
      
      # Get shared enrichments (if any)
      shared_enrichment_dataset = GSEA_enrichment_dataset %>%
                                  inner_join(ORA_enrichment_dataset, by = 'ID')
      
      module_enrichment[[dataset]] = shared_enrichment_dataset
    }
    
    top_modules_enrichment[[module]] = module_enrichment
  }
  
  return(top_modules_enrichment)
  
  
}

classification_model_dataset = function(dataset, modules_dataset, running_on_drake){
  
  # Extract information from input
  datExpr = dataset[['datExpr']]
  datMeta = dataset[['datMeta']]
  genes_info = dataset[['genes_info']]
  
  # Calculate Module-Diagnosis Correlation
  module_diagnosis_cor = module_diagonsis_association(datExpr, datMeta, modules_dataset)
  
  # Calculate Gene Significance
  GS_info = gene_significance(datExpr, datMeta)
  
  # Calculate Module Membership
  MM = module_membership(datExpr, datMeta, modules_dataset, running_on_drake)
  
  # Join results, create SFARI column and filter genes that weren't assigned a Module (gray module)
  classification_dataset = genes_info %>% 
    mutate(SFARI = !`gene-score` %in% c('Neuronal','Others')) %>% 
    dplyr::select(ID, SFARI) %>%
    left_join(modules_dataset, by = 'ID') %>%
    left_join(module_diagnosis_cor, by = 'Module') %>%
    left_join(GS_info, by='ID') %>%
    left_join(MM, by='ID') %>%
    filter(Module != 'gray') %>%
    dplyr::select(-Module) %>%
    relocate(SFARI, .after = last_col())
  
  # Move ID information to rownames
  rownames_dataset = classification_dataset$ID
  classification_dataset = classification_dataset %>% dplyr::select(-ID)
  rownames(classification_dataset) = rownames_dataset
  
  return(classification_dataset)
  
}
  
classification_model = function(Gandal_dataset, classification_dataset, SFARI_dataset, seed, correct_bias){
  
  # Extract expression matrix from input
  datExpr = Gandal_dataset[['datExpr']]
  datGenes = Gandal_dataset[['datGenes']]
  
  # CALCULATE PARAMETER TO CORRECT BIAS
  
  if(correct_bias){
    
    # Parameters
    p = 0.75
    Loops = 30
    
    # Run model to learn the lambda parameter
    model_output = run_weights_model(datExpr, classification_dataset, SFARI_dataset, p, seed, Loops)
    
    # Optimised bias correction parameter
    lambda = model_output[['lambda']]
    
  } else { lambda = 0 }
  
  
  # RUN FINAL MODEL
  
  # Parameters
  p = 0.75
  n_iter = 100
  seeds = seed:(seed+n_iter-1)
  
  # Run final model
  final_model_output = run_final_model_wrapper(datExpr, datGenes, classification_dataset, SFARI_dataset, 
                                               p, seeds, n_iter, lambda)
  
  if(correct_bias){# Add information from the bias correction process
    
    final_model_output = c(final_model_output, 'lambda' = lambda, 'BiasVector' = list(model_output$bias_vec),
                           'BalancedAccuracyVector' = list(model_output$b_acc_vec))
    
  }
  
  
  return(final_model_output)

} 

load_transcriptomic_dataset = function(dataset_name){

  dataset_path = 'Results/'
  if(dataset_name != 'Gandal') dataset_path = paste0('supportingDatasets/', dataset_name, '/' , dataset_path)
  
  load(paste0(dataset_path,'preprocessed_data.RData'))
  classification_dataset = read.csv(paste0(dataset_path,'classification_dataset.csv'), row.names=1)
  modules_dataset = read.csv(paste0(dataset_path,'modules_dataset.csv'))
  load(paste0(dataset_path,'biased_classification_model.RData'))
  load(paste0(dataset_path,'unbiased_classification_model.RData'))
  
  return(list('preprocessed_data'=preprocessed_data, 'classification_dataset'=classification_dataset,
              'modules_dataset'=modules_dataset,'biased_classification_model'=biased_classification_model,
              'unbiased_classification_model'=unbiased_classification_model))
}
