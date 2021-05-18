
# AUXILIARY FUNCTIONS

################################################################################################################
# PREPROCESSING
gg_colour_hue = function(n) {
  hues = seq(15, 375, length = n+1)
  pal = hcl(h = hues, l = 65, c = 100)[1:n]
}

SFARI_colour_hue = function(r) {
  pal = c('#FF7631','#FFB100','#E8E328','#8CC83F','#62CCA6','#59B9C9','#b3b3b3','#808080','gray','#d9d9d9')[r]
}

annotate_genes = function(datExpr, NCBI_dataset){
  
  # 1. Query archive version
  getinfo = c('ensembl_gene_id','external_gene_id','chromosome_name','start_position','end_position','strand')
  mart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', 
                 host = 'feb2014.archive.ensembl.org')
  datGenes = getBM(attributes = getinfo, filters=c('ensembl_gene_id'), values=rownames(datExpr), mart=mart) %>% 
    dplyr::rename('hgnc_symbol' = external_gene_id) %>%
    mutate(length = end_position-start_position)
  
  
  # 2. Get Biotype Labels
  # 2.1 Add NCBI annotations
  datGenes = datGenes %>% left_join(NCBI_dataset, by=c('ensembl_gene_id','hgnc_symbol'))
  # 2.2 Query current BioMart version for gene_biotype using Ensembl ID as key
  getinfo = c('ensembl_gene_id','gene_biotype')
  mart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl',
                 host = 'jan2020.archive.ensembl.org')
  datGenes_biotype = getBM(attributes = getinfo, filters = c('ensembl_gene_id'), mart=mart, 
                           values = datGenes$ensembl_gene_id[is.na(datGenes$gene_biotype)])
  datGenes = datGenes %>% 
    left_join(datGenes_biotype, by='ensembl_gene_id') %>%
    mutate(gene_biotype = coalesce(as.character(gene_biotype.x), gene_biotype.y)) %>%
    dplyr::select(-gene_biotype.x, -gene_biotype.y)
  
  
  # 3. Query current BioMart version for gene_biotype using gene symbol as key
  missing_genes = unique(datGenes$hgnc_symbol[is.na(datGenes$gene_biotype)])
  getinfo = c('hgnc_symbol','gene_biotype')
  datGenes_biotype_by_gene = getBM(attributes=getinfo, filters=c('hgnc_symbol'), mart=mart,
                                   values=missing_genes)
  dups = unique(datGenes_biotype_by_gene$hgnc_symbol[duplicated(datGenes_biotype_by_gene$hgnc_symbol)])
  datGenes_biotype_by_gene = datGenes_biotype_by_gene %>% filter(!hgnc_symbol %in% dups)
  datGenes = datGenes %>% 
    left_join(datGenes_biotype_by_gene, by='hgnc_symbol') %>% 
    mutate(gene_biotype = coalesce(gene_biotype.x, gene_biotype.y)) %>%
    dplyr::select(-gene_biotype.x, -gene_biotype.y)  
  
  
  # 4. Query feb2014 BioMart version for the missing biotypes
  missing_ensembl_ids = unique(datGenes$ensembl_gene_id[is.na(datGenes$gene_biotype)])
  getinfo = c('ensembl_gene_id','gene_biotype')
  mart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', 
                 host = 'feb2014.archive.ensembl.org')
  datGenes_biotype_archive = getBM(attributes = getinfo, filters=c('ensembl_gene_id'), 
                                   values = missing_ensembl_ids, mart=mart)
  datGenes = datGenes %>% 
    left_join(datGenes_biotype_archive, by='ensembl_gene_id') %>% 
    mutate(gene_biotype = coalesce(gene_biotype.x, gene_biotype.y)) %>%
    dplyr::select(-gene_biotype.x, -gene_biotype.y)
  
  # Reorder rows to match datExpr
  datGenes = datGenes[match(rownames(datExpr), datGenes$ensembl_gene_id),]
  
  return(datGenes)
  
}

filter_genes_and_samples = function(datExpr, datGenes, datMeta, threshold){
  
  # 1. Filter entries that don't correspond to genes
  to_keep = !is.na(datGenes$length)
  datGenes = datGenes[to_keep,]
  datExpr = datExpr[to_keep,]
  rownames(datGenes) = rownames(datExpr)
  
  # 2. Filter genes that do not encode any protein
  to_keep = datGenes$gene_biotype=='protein_coding'
  datGenes = datGenes[to_keep,]
  datExpr = datExpr[to_keep,]
  
  # 3. Filter genes with low expression levels (using high percentage of zeros to measure this)
  # threhsold = Minimum percentage of non-zero entries allowed per gene
  to_keep = apply(datExpr, 1, function(x) 100*mean(x>0)) >= threshold
  datGenes = datGenes[to_keep,]
  datExpr = datExpr[to_keep,]
  
  # 4. Filter outlier samples
  # 4.1 Remove samples belonging to Subject AN03345 (Following Gandal's script)
  to_keep = (datMeta$Subject_ID != 'AN03345')
  datMeta = datMeta[to_keep,]
  datExpr = datExpr[,to_keep]
  # 4.2 Remove samples farther away than 2SD from the rest
  absadj = datExpr %>% bicor %>% abs
  netsummary = fundamentalNetworkConcepts(absadj)
  ku = netsummary$Connectivity
  z.ku = (ku-mean(ku))/sqrt(var(ku))
  to_keep = z.ku > -2
  datMeta = datMeta[to_keep,]
  datExpr = datExpr[,to_keep]
  
  # 5. Filter repeated genes
  dup_genes = datGenes$hgnc_symbol %>% duplicated
  datGenes = datGenes[!dup_genes,]
  datExpr = datExpr[!dup_genes,]
  
  filtered_data = list('datExpr'=datExpr, 'datGenes'=datGenes, 'datMeta'=datMeta)
  
  return(filtered_data)
}

perform_sva = function(dataset){
  
  datExpr = dataset[['datExpr']] %>% data.frame
  datGenes = dataset[['datGenes']] %>% data.frame
  datMeta = dataset[['datMeta']]
  rm(dataset)
  
  # Create a DeseqDataSet object, estimate the library size correction and save the normalized counts matrix
  counts = datExpr %>% as.matrix
  rowRanges = GRanges(datGenes$chromosome_name,
                      IRanges(datGenes$start_position, width=datGenes$length),
                      strand=datGenes$strand,
                      feature_id=datGenes$ensembl_gene_id)
  se = SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=rowRanges, colData=datMeta)
  dds = DESeqDataSet(se, design = ~ Diagnosis)
  
  dds = estimateSizeFactors(dds)
  norm.cts = counts(dds, normalized = TRUE)
  
  # Provide the normalized counts and two model matrices to SVA
  mod = model.matrix(~ Diagnosis, colData(dds))
  mod0 = model.matrix(~ 1, colData(dds))
  sva_fit = svaseq(norm.cts, mod=mod, mod0=mod0)
  
  # Include SV estimations to datMeta information
  sv_data = sva_fit$sv %>% data.frame
  colnames(sv_data) = paste0('SV', 1:ncol(sv_data))
  datMeta = cbind(datMeta, sv_data)
  
  return(datMeta)
  
}

normalisation_and_DEA = function(datExpr, datGenes, datMeta){
  
  # Create DESeq dataset including batch variables and Diagnosis into the design
  counts = datExpr %>% as.matrix
  rowRanges = GRanges(datGenes$chromosome_name,
                      IRanges(datGenes$start_position, width=datGenes$length),
                      strand=datGenes$strand,
                      feature_id=datGenes$ensembl_gene_id)
  se = SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=rowRanges, colData=datMeta)
  dds = DESeqDataSet(se, design = ~ Batch + SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + 
                       SV10 + SV11 + SV12 + SV13 + Diagnosis)
  
  # Perform DEA
  dds = DESeq(dds)
  DE_info = results(dds, lfcThreshold=0, altHypothesis='greaterAbs')
  
  # Perform normalisation using vst
  vsd = vst(dds)
  
  # Extract data
  datExpr = assay(vsd)
  datMeta = colData(vsd) %>% data.frame
  
  norm_data = list('datExpr'=datExpr, 'datMeta'=datMeta, 'DE_info'=DE_info, 'dds'=dds)
  
  return(norm_data)
}

correct_batch_effects = function(dataset){
  
  datExpr = dataset[['datExpr']] %>% data.frame
  datMeta = dataset[['datMeta']]
  dds = dataset[['dds']]
  rm(dataset)
  
  
  # SVA
  correctDatExpr = function(datExpr, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(datExpr))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(datExpr - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
  }
  
  mod = model.matrix(~ Diagnosis, colData(dds))
  svs = datMeta %>% dplyr::select(SV1:SV13) %>% as.matrix
  datExpr = correctDatExpr(as.matrix(datExpr), mod, svs)
  
  # ComBat
  datExpr = datExpr %>% as.matrix %>% ComBat(batch=datMeta$Batch)
  
  return(datExpr)
}

################################################################################################################
# WGCNA
module_diagonsis_association = function(datExpr, datMeta, modules_dataset){
  
  # Extract Diagnosis information from MetaData
  datTraits = datMeta %>% dplyr::select(Diagnosis)
  
  # Recalculate MEs with color labels
  ME_object = datExpr %>% t %>% moduleEigengenes(colors = modules_dataset$Module)
  MEs = orderMEs(ME_object$eigengenes)
  
  # Calculate correlation between eigengenes and Diagnosis and their p-values
  moduleTraitCor = MEs %>% apply(2, function(x) hetcor(x, datTraits)$correlations[1,-1])
  module_diagnosis_cor = data.frame('Module' = gsub('ME','',colnames(MEs)),
                                    'MTcor' = moduleTraitCor)
  
  return(module_diagnosis_cor)
}

gene_significance = function(datExpr, datMeta){
  
  # Calculate Gene Significance (correlation between gene expression profile and Diagnosis)
  GS_info = data.frame('ID'=rownames(datExpr),
                       'GS'=datExpr %>% apply(1, function(x) hetcor(x,datMeta$Diagnosis)$correlations[1,2]))
  
  # Complete missing Gene Significance values in dataset
  GS_missing = GS_info$ID[is.na(GS_info$GS)] %>% as.character
  if(length(GS_missing)>0){
    for(g in GS_missing) GS_info$GS[GS_info$ID == g] = polyserial(as.numeric(datExpr[g,]), datMeta$Diagnosis)
  }
  
  # Add extra variable with the magnitude of the GS
  GS_info = GS_info %>% mutate(absGS = abs(GS))
  
  return(GS_info)
}

module_membership = function(datExpr, datMeta, modules_dataset, running_on_drake){
  
  # Check if the MM information is already stored in IntermediateData/ to avoid running it again
  run_MM = TRUE
  mm_file_name = './IntermediateData/MM.csv'
  
  if(file.exists(mm_file_name)){
    MM = read.csv(mm_file_name)
    modules = gsub('#','MM.', modules_dataset$Module[modules_dataset$Module!='gray'] %>% unique)
    if(nrow(MM) == nrow(datExpr) & ncol(MM) == length(modules)+1 & all(modules %in% as.character(colnames(MM)))) {
      run_MM = FALSE
    }
  }
  
  # This can take a while to run, it's preferable to run in a computer with many cores
  if(run_MM){
    
    # Extract Diagnosis information from MetaData
    datTraits = datMeta %>% dplyr::select(Diagnosis)
    
    # Recalculate MEs with color labels
    ME_object = datExpr %>% t %>% moduleEigengenes(colors = modules_dataset$Module)
    MEs = orderMEs(ME_object$eigengenes)
    
    if(running_on_drake){
      
      MM = apply(MEs, 2, function(x) hetcor(as.numeric(datExpr[i,]), x)$correlations[1,2])
      
    } else {
      
      #########################################################################################################
      # NOTE:
      # This implementation is much faster than the one if running_on_drake = TRUE but drake produces errors
      # when used together with the parallel library https://github.com/ropensci/drake/issues/925
      #########################################################################################################
      
      # Load doParallel library
      library(doParallel)
      
      #setup parallel backend to use many processors
      cores = detectCores()
      cl = makeCluster(cores-1)
      registerDoParallel(cl)
      
      # Create matrix with MM by gene
      MM = foreach(i=1:nrow(datExpr), .combine=rbind) %dopar% {
        library(polycor)
        tempMatrix = apply(MEs, 2, function(x) hetcor(as.numeric(datExpr[i,]), x)$correlations[1,2])
        tempMatrix
      }
      
      # Stop clusters
      stopCluster(cl)
      
    }

    rownames(MM) = rownames(datExpr)
    colnames(MM) = paste0('MM',gsub('ME','',colnames(MEs)))
    MM = MM %>% as.data.frame %>% mutate(ID = rownames(.)) %>% dplyr::select(-MMgray)
    
    # Save MM
    write.csv(MM, file = mm_file_name, row.names = FALSE)
  }

  return(MM)
  
}

prepare_dataset_for_EA = function(modules_dataset){
  
  EA_dataset = modules_dataset %>% 
               dplyr::rename('ensembl_gene_id' = ID) %>% 
               filter(Module!='gray')
  
  # ClusterProfile works with Entrez Gene Ids, o we have to assign one to each gene
  getinfo = c('ensembl_gene_id','entrezgene')
  mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', 
                 host='feb2014.archive.ensembl.org')
  biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), 
                         values=EA_dataset$ensembl_gene_id, mart=mart)
  
  EA_dataset = biomart_output %>% 
               left_join(EA_dataset, by = 'ensembl_gene_id') %>% 
               dplyr::rename('ID' = ensembl_gene_id)
               
  
  return(EA_dataset)
  
}

enrichment_with_GSEA = function(EA_dataset, top_modules, classification_dataset, nPerm = 1e5){
  
  EA_dataset = EA_dataset %>% dplyr::select(-Module) %>%
               left_join(classification_dataset %>% mutate(ID = rownames(.)) %>% 
                         dplyr::select(ID, contains('MM.')), by='ID')
  
  GSEA_enrichment = list()
  
  for(module in top_modules){
    
    cat(paste0('\nModule: ', which(top_modules == module), '/', length(top_modules)))
    
    geneList = EA_dataset %>% pull(paste0('MM.',substring(module,2)))
    names(geneList) = EA_dataset %>% pull(entrezgene) %>% as.character
    geneList = sort(geneList, decreasing = TRUE)
    
    GSEA_GO = gseGO(geneList, OrgDb = org.Hs.eg.db, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                    nPerm = nPerm, verbose = FALSE, seed = TRUE)
    
    GSEA_DO = gseDO(geneList, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                    nPerm = nPerm, verbose = FALSE, seed = TRUE)
    
    GSEA_DGN = gseDGN(geneList, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                      nPerm = nPerm, verbose = FALSE, seed = TRUE)
    
    GSEA_KEGG = gseKEGG(geneList, organism = 'human', pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                        nPerm = nPerm, verbose = FALSE, seed = TRUE)
    
    GSEA_Reactome = gsePathway(geneList, organism = 'human', pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                               nPerm = nPerm, verbose = FALSE, seed = TRUE)
    
    GSEA_enrichment[[module]] = list('GO' = GSEA_GO, 'DO' = GSEA_DO, 'DGN' = GSEA_DGN, 'KEGG' = GSEA_KEGG, 'Reactome' = GSEA_Reactome)
  }
  
  return(GSEA_enrichment)
  
}

enrichment_with_ORA = function(EA_dataset, top_modules){
  
  # Prepare input
  universe = EA_dataset$entrezgene %>% as.character
  
  # Perform Enrichment
  
  ORA_enrichment = list()
  
  for(module in top_modules){
    
    genes_in_module = EA_dataset %>% filter(Module == module) %>% pull(entrezgene)
    
    ORA_GO = enrichGO(gene = genes_in_module, universe = universe, OrgDb = org.Hs.eg.db, ont = 'All', 
                      pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, qvalueCutoff = 1)
    
    ORA_DO = enrichDO(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                      pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_DGN = enrichDGN(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                        pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_KEGG = enrichKEGG(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                          pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1) 
    
    ORA_Reactome = enrichPathway(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                                 pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_enrichment[[module]] = list('GO' = ORA_GO, 'DO' = ORA_DO, 'DGN' = ORA_DGN, 'KEGG' = ORA_KEGG, 'Reactome' = ORA_Reactome)
  }
  
  return(ORA_enrichment)
  
}

################################################################################################################
# CLASSIFICATION MODEL

create_train_test_sets = function(dataset, SFARI_dataset, p, seed){
  
  # Get SFARI Score of all the samples so our train and test sets are balanced for each score
  sample_scores = data.frame('ID' = rownames(dataset)) %>% 
                  left_join(SFARI_dataset %>% dplyr::select(ID,`gene-score`), by = 'ID') %>% 
                  mutate(gene.score = ifelse(is.na(`gene-score`), 'None', `gene-score`)) %>%
                  dplyr::select(-`gene-score`)
  
  set.seed(seed)
  train_idx = createDataPartition(sample_scores$gene.score, p = p, list = FALSE)
  train_set = dataset[train_idx,]
  test_set = dataset[-train_idx,]
  
  # Modify SFARI label in train set, save gene IDS (bc we lose them with SMOTE) and perform oversampling using SMOTE
  set.seed(seed)
  train_set = train_set %>% mutate(SFARI = ifelse(SFARI == TRUE, 'SFARI', 'not_SFARI') %>% as.factor,
                                   ID = rownames(.) %>% as.factor) %>% SMOTE(form = SFARI ~ . - ID)
  train_set_IDs = train_set %>% pull(ID)
  
  return(list('train_set' = train_set %>% dplyr::select(-ID), 'test_set' = test_set, 
              'train_set_IDs' = train_set_IDs))
}

calculate_performance_metrics = function(predictions){
  
  if(all(predictions$pred == 0)){
    prec = NA
    F1 = NA
  } else {
    prec = Precision(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
    F1 = F1_Score(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  }
  
  acc = mean(predictions$SFARI == predictions$pred)
  rec = Recall(predictions$SFARI %>% as.numeric, predictions$pred %>% as.numeric, positive = '1')
  pred_ROCR = prediction(predictions$prob, predictions$SFARI)
  AUC = performance(pred_ROCR, measure='auc')@y.values[[1]]
  MLP = performance(pred_ROCR, measure='lift', x.measure='rpp')@y.values[[1]] %>% max(na.rm = TRUE)
  b_acc = mean(c(mean(predictions$SFARI[predictions$SFARI] == predictions$pred[predictions$SFARI]),
                 mean(predictions$SFARI[!predictions$SFARI] == predictions$pred[!predictions$SFARI])))
  
  return(list('acc' = acc, 'prec' = prec, 'rec' = rec, 'F1' = F1, 'AUC' = AUC, 'MLP'=MLP, 'b_acc' = b_acc))
  
}

run_weights_model = function(datExpr, classification_dataset, SFARI_dataset, p, seed, Loops){
  
  # CREATE TRAIN AND TEST SETS
  train_test_sets = create_train_test_sets(classification_dataset, SFARI_dataset, p, seed)
  train_set = train_test_sets[['train_set']]
  test_set = train_test_sets[['test_set']]
  train_set_IDs = train_test_sets[['train_set_IDs']]
  
  
  # SET INITIAL PARAMETERS
  
  # General parameters
  lambda_seq = 10^seq(1, -4, by = -.1)
  k_fold = 10
  cv_repeats = 5
  set.seed(seed)
  trControl = trainControl(method = 'repeatedcv', number = k_fold, repeats = cv_repeats, verboseIter = FALSE, 
                           classProbs = TRUE, savePredictions = 'final', summaryFunction = twoClassSummary,
                           seeds = as.list(seq(seed, seed+length(lambda_seq))))
  # Bias correction parameters
  eta = 0.5
  lambda = 0
  w = rep(1, nrow(train_set))
  
  
  # TRAIN MODEL
  set.seed(seed)
  h = train(SFARI ~., data = train_set, method = 'glmnet', trControl = trControl, metric = 'ROC',
            tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  
  
  # CORRECT BIAS
  
  # Mean Expression info
  mean_expr = data.frame('ID' = train_set_IDs) %>% 
              left_join(data.frame('ID' = rownames(datExpr), 'meanExpr' = rowMeans(datExpr)), by = 'ID') %>%
              mutate('meanExpr_std' = (meanExpr-mean(meanExpr))/sd(meanExpr))
  
  # Track behaviour of plot
  bias_vec = c()
  b_acc_vec = c()
  
  for(l in 1:Loops){
    
    # Calculate bias for positive predicted samples
    bias = mean(mean_expr$meanExpr_std[predict(h, train_set)=='SFARI'])
    if(is.na(bias)) bias = 0 # This happens when all the observations are labelled Negative
    
    # Update weights
    lambda = lambda - eta*bias
    w_hat = exp(lambda*mean_expr$meanExpr_std)
    w = 1/(1+w_hat)
    w[train_set$SFARI=='SFARI'] = w[train_set$SFARI=='SFARI']*w_hat[train_set$SFARI=='SFARI']
    
    # Update tracking vars
    bias_vec = c(bias_vec, bias)
    preds = predict(h, train_set)
    b_acc = mean(c(mean(preds[preds=='SFARI'] == train_set$SFARI[preds=='SFARI']),
                   mean(preds[preds!='SFARI'] == train_set$SFARI[preds!='SFARI'])))
    b_acc_vec = c(b_acc_vec, b_acc)
    #acc_vec = c(acc_vec, mean(predict(h, train_set) == train_set$SFARI))
    
    # Update h
    set.seed(seed)
    h = train(SFARI ~., data = train_set, method = 'glmnet', weights = w, trControl = trControl, 
              metric = 'ROC', tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  }
  
  
  return(list('lambda' = lambda, 'bias_vec' = bias_vec, 'b_acc_vec' = b_acc_vec))
}

run_final_model = function(datExpr, classification_dataset, SFARI_dataset, p, seed, lambda){
  
  # CREATE TRAIN AND TEST SETS
  train_test_sets = create_train_test_sets(classification_dataset, SFARI_dataset, p, seed)
  train_set = train_test_sets[['train_set']]
  test_set = train_test_sets[['test_set']]
  train_set_IDs = train_test_sets[['train_set_IDs']]
  
  
  # SET INITIAL PARAMETERS
  
  # General parameters
  lambda_seq = 10^seq(1, -4, by = -.1)
  k_fold = 10
  cv_repeats = 5
  set.seed(seed)
  trControl = trainControl(method = 'repeatedcv', number = k_fold, repeats = cv_repeats, verboseIter = FALSE, 
                           classProbs = TRUE, savePredictions = 'final', summaryFunction = twoClassSummary,
                           seeds = as.list(seq(seed*100, seed*100+length(lambda_seq))))
  
  # Bias correcting parameters
  mean_expr = data.frame('ID' = train_set_IDs) %>% 
              left_join(data.frame('ID' = rownames(datExpr), 'meanExpr' = rowMeans(datExpr)), by = 'ID') %>%
              mutate('meanExpr_std' = (meanExpr-mean(meanExpr))/sd(meanExpr))
  w_hat = exp(lambda*mean_expr$meanExpr_std)
  w = 1/(1+w_hat)
  w[train_set$SFARI=='SFARI'] = w[train_set$SFARI=='SFARI']*w_hat[train_set$SFARI=='SFARI']
  
  
  # TRAIN MODEL
  set.seed(seed)
  fit = train(SFARI ~., data = train_set, method = 'glmnet', weights = w, trControl = trControl, 
              metric = 'ROC', tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq))
  
  
  # PREDICT TEST SET LABELS AND CREATE PERFORMANCE METRICS
  
  # Predict labels in test set
  predictions = fit %>% predict(test_set, type = 'prob')
  preds = data.frame('ID' = rownames(test_set), 'prob' = predictions$SFARI) %>% 
          mutate(pred = prob > 0.5)
  
  
  # Measure performance of the model
  pm_dataset = preds %>% left_join(data.frame('ID' = rownames(test_set), 'SFARI' = test_set$SFARI), by = 'ID')
  performance_metrics = calculate_performance_metrics(pm_dataset)
  
  # Extract coefficients from features
  coefs = coef(fit$finalModel, fit$bestTune$lambda) %>% as.vector
  
  
  return(list('performance_metrics' = performance_metrics, 'preds' = preds, 'coefs'= coefs))
  
}

run_final_model_wrapper = function(datExpr, datGenes, classification_dataset, SFARI_dataset, p, seeds, n_iter, lambda){
  
  # Store outputs
  acc = c() ; prec = c() ; rec = c() ; F1 = c() ; AUC = c() ; MLP = c() ; b_acc = c()
  predictions = data.frame('ID' = rownames(classification_dataset), 'SFARI' = classification_dataset$SFARI, 
                           'prob' = 0, 'pred' = 0, 'n' = 0)
  coefs = data.frame('var' = c('Intercept', colnames(classification_dataset)[-ncol(classification_dataset)]),
                     'coef' = 0)
  
  for(seed in seeds){
    
    # Run model
    model_output = run_final_model(datExpr, classification_dataset, SFARI_dataset, p, seed, lambda)
    
    # Update outputs
    acc = c(acc, model_output$performance_metrics$acc)
    prec = c(prec, model_output$performance_metrics$prec)
    rec = c(rec, model_output$performance_metrics$rec)
    F1 = c(F1, model_output$performance_metrics$F1)
    AUC = c(AUC, model_output$performance_metrics$AUC)
    MLP = c(MLP, model_output$performance_metrics$MLP)
    b_acc = c(b_acc, model_output$performance_metrics$b_acc)
    preds = model_output$performance_metrics$preds
    coefs$coef = coefs$coef + model_output$coefs
    update_preds = model_output$preds %>% dplyr::select(-ID) %>% mutate(n=1)
    predictions[predictions$ID %in% model_output$preds$ID, c('prob','pred','n')] = 
      predictions[predictions$ID %in% model_output$preds$ID, c('prob','pred','n')] + update_preds
  }
  
  # Process results
  coefs = coefs %>% mutate(coef = coef/n_iter)
  
  predictions = predictions %>% mutate(prob = prob/n, pred_count = pred, pred = prob>0.5) %>%
                left_join(datGenes %>% mutate(ID=ensembl_gene_id) %>% dplyr::select(ID, hgnc_symbol), by = 'ID')
  
  pm_interim_models = list('acc' = acc, 'prec' = prec, 'rec' = rec, 'F1' = F1, 'AUC' = AUC, 'MLP' = MLP, 'b_acc' = b_acc)
  
  # Calculate performance metrics of the final model
  final_predictions = predictions %>% drop_na(pred) #filter(!is.na(pred))
  pm_final_model = calculate_performance_metrics(final_predictions)
  
  
  return(list('predictions' = predictions, 'coefficients' = coefs, 'pm_interim_models' = pm_interim_models, 
              'pm_final_model' = pm_final_model))
}

