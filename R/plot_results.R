
# PLOT RESULTS

source('R/packages.R')

SFARI_colour_hue = function(r) {
  pal = c('#FF7631','#FFB100','#E8E328','#808080','#b3b3b3')[r]
}

################################################################################################################
################################################################################################################
# LOAD DATA

# SFARI Genes
SFARI_dataset = read.csv('Results/SFARI_dataset.csv') %>% dplyr::rename('gene-score' = gene.score)

# Transcriptomic data
dataset_name = 'Gandal'
transcriptomic_data = load_transcriptomic_dataset(dataset_name)
preprocessed_data = transcriptomic_data$preprocessed_data
datExpr = preprocessed_data$datExpr
datMeta = preprocessed_data$datMeta
genes_info = preprocessed_data$genes_info
dds = preprocessed_data$dds # Wright's preprocessed data doesn't contain this last element
classification_dataset = transcriptomic_data$classification_dataset
modules_dataset = transcriptomic_data$modules_dataset
biased_classification_model = transcriptomic_data$biased_classification_model
unbiased_classification_model = transcriptomic_data$unbiased_classification_model

# Performance metrics
dataset_path = 'Results/'
if(dataset_name != 'Gandal') dataset_path = paste0('supportingDatasets/', dataset_name, '/' , dataset_path)
pm_unbiased_model = read.csv(paste0(dataset_path,'performance_metrics_unbiased_model.csv'))
pm_biased_model = read.csv(paste0(dataset_path,'performance_metrics_biased_model.csv'))
pm_shuffled_model = read.csv(paste0(dataset_path,'performance_metrics_shuffled_labels_model.csv'))

# Other scoring systems
DisGeNET = list()
DisGeNET[['asd']] = disease2gene(disease = 'C0004352')@qresult   # Autism Spectrum Disorder
DisGeNET[['scz']] = disease2gene(disease = 'C0036341')@qresult   # Schizophrenia
DisGeNET[['bd']]  = disease2gene(disease = 'C0005586')@qresult   # Bipolar Disorder
DisGeNET[['id']]  = disease2gene(disease = 'C3714756')@qresult   # Intellectual Disability
DisGeNET[['dd']]  = disease2gene(disease = 'C0011581')@qresult   # Depressive Disorder
DisGeNET[['ai']]  = disease2gene(disease = 'C0001973')@qresult   # Alcoholic Intoxication, Chronic

krishnan = read_excel('InputData/krishnan_probability_score.xlsx', sheet = 'Genome-wide_prediction') %>%
           dplyr::select(symbol, score, `p-value`) %>% 
           dplyr::rename('gene' = symbol, 'Krishnan' = score, 'Krishnan_p_value' = `p-value`)

TADA = read_excel('InputData/sanders_TADA_score.xlsx') %>% dplyr::select(Gene, p.TADA, q.TADA) %>%
       dplyr::rename('gene' = Gene, 'TADA' = q.TADA, 'TADA_p_value' = p.TADA)


################################################################################################################
################################################################################################################
# PLOTS FOR THE RESULTS SECTION

################################################################################################################
# SECTION 1. SFARI genes have higher levels of expression than other neuronal and non-neuronal genes

# Fig 1. Mean level of expression plays a central role in gene characterisation
# Size: 7 x 6 in and then crop for main figure, 3 x 7 for supplementary figures
pca = datExpr %>% prcomp
pca_genes = data.frame('PC1' = pca$x[,1], 'PC2' = pca$x[,2], 'MeanExpression' = rowMeans(datExpr)) %>% #sample_n(10000) %>%
            ggplot(aes(PC1, PC2, color=MeanExpression)) + geom_point(alpha=0.3) + theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom') + coord_fixed() +
            scale_color_viridis() + scale_y_continuous(breaks=c(-4,0,4), limits = c(-5,5)) + 
            guides(color = guide_colourbar(barwidth = 10, barheight = 0.5, ticks = F, title = 'Mean Expression ')) +
            xlab(paste0('PC1 (',round(100*summary(pca)$importance[2,1],1),'%)')) +
            ylab(paste0('PC2 (',round(100*summary(pca)$importance[2,2],1),'%)'))

#ggarrange(pca_genes_Gupta, pca_genes_Wright, nrow=2, common.legend = TRUE, legend='bottom')

# Fig 2. SFARI genes have higher levels of expression than other genes
# Size: 3.5 x 6 in form main figure, 7 x 7 for supplementary figures

# Fig 2A. SFARI Genes
plot_data = data.frame('ID'=rownames(datExpr), 'MeanExpr'=rowMeans(datExpr)) %>% 
            left_join(genes_info, by='ID') %>% 
            mutate(Group = factor(ifelse(gene.score %in% c('Neuronal','Others'), gene.score, 'SFARI'), 
                   levels = c('SFARI', 'Neuronal', 'Others')))
wt = plot_data %>% wilcox_test(MeanExpr~Group, p.adjust.method='BH') %>% add_x_position(x = 'group')
increase = 1 ; base = 15.7
pos_y_comparisons = c(base, base+increase, base)
mean_expr_SFARI = plot_data %>% ggplot(aes(Group, MeanExpr)) + 
                  geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, aes(fill=Group)) + 
                  stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .03, coord.flip = TRUE) +
                  scale_fill_manual(values=c('#00A4F7', SFARI_colour_hue(r=4:5))) + xlab('') + 
                  ylab('Mean Expression') + theme_minimal() + 
                  scale_x_discrete(labels = c('SFARI\ngenes', 'Neuronal\ngenes', 'Other\ngenes')) +
                  theme(legend.position='none')#, plot.margin = unit(c(1,200,1,1), 'points'))

# Fig 2B. SFARI Scores
plot_data = data.frame('ID'=rownames(datExpr), 'MeanExpr'=rowMeans(datExpr)) %>% 
            left_join(genes_info, by='ID')
wt = plot_data %>% wilcox_test(MeanExpr~gene.score, p.adjust.method='BH') %>% add_x_position(x = 'group')
increase = 1 ; base = 15.5
pos_y_comparisons = c(base + c(0,1,3,5)*increase, base+c(0,2,4)*increase, base+0:1*increase, base)
mean_expr_SFARI_scores = plot_data %>% ggplot(aes(gene.score, MeanExpr)) + 
                          geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, 
                                       aes(fill=gene.score)) +
                          stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .01) +
                          scale_fill_manual(values=SFARI_colour_hue(r=1:5)) + xlab('') + ylab('Mean Expression') + 
                          scale_x_discrete(labels = c('SFARI\nScore 1', 'SFARI\nScore 2', 'SFARI\nScore 3',
                                                      'Neuronal\ngenes','Other\ngenes')) +
                          theme_minimal() + theme(legend.position='none')#, plot.margin = unit(c(1,200,1,1), 'points'))

grid.arrange(mean_expr_SFARI_Gupta, mean_expr_SFARI_scores_Gupta, mean_expr_SFARI_Wright, mean_expr_SFARI_scores_Wright, 
             nrow = 2, widths = c(0.4,0.6))


rm(plot_data, wt, increase, base, pos_y_comparisons)

################################################################################################################
# SECTION 2. Gene level: SFARI genes have smaller differences in level of expression between ASD and control 
# patients than other neuronal genes

# Fig 3. SFARI Genes have a lower percentage of differentially expressed genes than neuronal genes and a similar 
# percentage to the rest of the genes
# Size: 4 x 5 in for main figure and for supplementary figures
lfc_list = seq(1, 1.2, 0.01)
genes_info = genes_info %>% mutate(Group=factor(ifelse(gene.score %in% c('Neuronal','Others'), gene.score, 'SFARI'), 
                                                levels = c('SFARI', 'Neuronal', 'Others')))
all_counts = data.frame('group'='All', 'n'=as.character(nrow(genes_info)))
Others_counts = data.frame('group'='Others', n=as.character(sum(genes_info$Group == 'Others')))
Neuronal_counts = data.frame('group'='Neuronal', n=as.character(sum(genes_info$Neuronal)))
lfc_counts_all = genes_info %>% filter(Group=='SFARI') %>% tally %>% mutate('group'='SFARI', 'n'=as.character(n)) %>%
                 dplyr::select(group, n) %>% bind_rows(Neuronal_counts, Others_counts, all_counts) %>%
                 mutate('lfc'=-1) %>%  dplyr::select(lfc, group, n)

for(lfc in lfc_list){
  # Recalculate genes_info with the new threshold (p-values change)
  DE_genes = results(dds, lfcThreshold=log2(lfc), altHypothesis='greaterAbs') %>% data.frame %>% mutate('ID'=rownames(.)) %>% 
             left_join(genes_info %>% dplyr::select(ID, Neuronal, gene.score, Group), by = 'ID') %>% 
             filter(padj<0.05 & abs(log2FoldChange)>log2(lfc))
  # Calculate counts by groups
  all_counts = data.frame('group'='All', 'n'=as.character(nrow(DE_genes)))
  Others_counts = data.frame('group'='Others', n=as.character(sum(DE_genes$Group == 'Others')))
  Neuronal_counts = data.frame('group'='Neuronal', n=as.character(sum(DE_genes$Neuronal)))
  lfc_counts = DE_genes %>% filter(Group == 'SFARI') %>% tally %>% mutate('group'='SFARI', 'n'=as.character(n)) %>%
               bind_rows(Neuronal_counts, Others_counts, all_counts) %>% mutate('lfc'=lfc) %>% dplyr::select(lfc, group, n)
  # Update lfc_counts_all
  lfc_counts_all = lfc_counts_all %>% bind_rows(lfc_counts)
}
# Add missing entries with 0s
lfc_counts_all = expand.grid('group'=unique(lfc_counts_all$group), 'lfc'=unique(lfc_counts_all$lfc)) %>% 
                 left_join(lfc_counts_all, by=c('group','lfc')) %>% replace(is.na(.), 0)
# Calculate percentage of each group remaining
tot_counts = genes_info %>% filter(Group == 'SFARI') %>% tally() %>% mutate('group'='SFARI', 'tot'=n) %>% 
             dplyr::select(group, tot) %>%
             bind_rows(data.frame('group'='Neuronal', 'tot'=sum(genes_info$Neuronal)),
                       data.frame('group' = 'Others', 'tot' = sum(genes_info$Group == 'Others')),
                       data.frame('group'='All', 'tot'=nrow(genes_info)))
lfc_counts_all = lfc_counts_all %>% filter(lfc!=-1) %>% #, group!='Others') %>% 
                 left_join(tot_counts, by='group') %>% mutate('perc'=round(100*as.numeric(n)/tot,2))

perc_DE_genes = lfc_counts_all %>% filter(group != 'All') %>% 
                mutate(group = factor(group, levels = c('SFARI', 'Neuronal', 'Others'))) %>%
                ggplot(aes(log(lfc), perc, color = group, shape = group)) + geom_point() + geom_line() +
                scale_color_manual(values=c('#00A4F7', SFARI_colour_hue(r=c(4,5)))) + 
                ylab('% of differentially expressed genes') +  xlab('Log Fold Change') + 
                labs(color = 'Group', shape = 'Group') + theme_minimal() + theme(legend.position = 'bottom')


# Fig 4. SFARI genes have lower log fold-change magnitudes than neuronal genes and similar magnitudes to non-neuronal genes
# Size: 3.5 x 6 in for main figure and 7 x 7 for supplementary figures

#DE_info = results(dds, lfcThreshold=0, altHypothesis='greaterAbs')
#DE_info_shrunken = lfcShrink(dds, coef='Diagnosis_ASD_vs_CTL', lfcThreshold=0, res = DE_info, quiet=TRUE)
#shrunken_lfc_info = DE_info_shrunken %>% data.frame %>% mutate(ID = rownames(DE_info_shrunken %>% data.frame))

# Fig 4A. SFARI Genes
plot_data = data.frame('ID'=rownames(datExpr)) %>% #left_join(shrunken_lfc_info, by = 'ID') %>%
            left_join(genes_info %>% dplyr::select(ID, gene.score, log2FoldChange, shrunken_log2FoldChange), by='ID') %>%
            mutate(Group = factor(ifelse(gene.score %in% c('Neuronal','Others'), gene.score, 'SFARI'), 
                                  levels = c('SFARI', 'Others', 'Neuronal')))
wt = plot_data %>% mutate(abs_lfc = abs(shrunken_log2FoldChange)) %>% 
     wilcox_test(abs_lfc~Group, p.adjust.method='BH') %>% add_x_position(x = 'group')
increase = 0.03 # 0.006 for Wright
base = 0.35 # 0.075 for Wright
pos_y_comparisons = c(base, base + increase, base)
lfc_magnitude = plot_data %>% mutate(Group = factor(Group, levels = c('SFARI','Others', 'Neuronal'))) %>% 
                ggplot(aes(Group, abs(shrunken_log2FoldChange))) + 
                geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, aes(fill=Group)) +
                stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .01) +
                scale_fill_manual(values=c('#00A4F7', SFARI_colour_hue(r=c(5,4)))) + 
                coord_cartesian(ylim= c(0, max(pos_y_comparisons))) +
                scale_x_discrete(labels = c('SFARI\ngenes', 'Other\ngenes', 'Neuronal\ngenes')) +
                xlab('') + ylab('Log fold-change magnitude') + theme_minimal() + 
                theme(legend.position='none')#, plot.margin = unit(c(1,200,1,1), 'points'))


# Fig 4B. SFARI Scores
wt = genes_info %>% mutate(abs_lfc = abs(shrunken_log2FoldChange), 
                           Group = factor(gene.score, levels = c('1','2','3','Others', 'Neuronal'))) %>% 
  wilcox_test(abs_lfc~Group, p.adjust.method='BH') %>% add_x_position(x = 'group')

increase = 0.03 # 0.01 (Wright)
base = 0.35 # 0.075 (Wright)
pos_y_comparisons = c(base + c(0,1,3,5)*increase, base+c(0,2,4)*increase, base+0:1*increase, base)
lfc_magnitude_by_score = genes_info %>% 
                         mutate(gene.score = factor(gene.score, levels = c('1','2','3','Others','Neuronal'))) %>%
                         ggplot(aes(gene.score, abs(shrunken_log2FoldChange))) + 
                         geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, aes(fill=gene.score)) + 
                         stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .01) + 
                         scale_fill_manual(values=SFARI_colour_hue(r=c(1:3,5,4))) + 
                         coord_cartesian(ylim= c(0, max(pos_y_comparisons))) + 
                         scale_x_discrete(labels=c('1' = 'SFARI\nScore 1', '2' = 'SFARI\nScore 2', '3' = 'SFARI\nScore 3', 
                                                  'Neuronal' = 'Neuronal\ngenes', 'Others' = 'Other\ngenes')) +
                         xlab('') + ylab('Log fold-change magnitude') + theme_minimal() + theme(legend.position='none')
                         #, plot.margin = unit(c(1,200,1,1), 'points'))

# grid.arrange(lfc_magnitude, lfc_magnitude_by_score, nrow = 1, widths = c(0.4, 0.6))


rm(lfc_magnitude, Neuronal_counts, Others_counts, lfc, lfc_list, plot_data, wt, increase, base, 
   pos_y_comparisons)

################################################################################################################
# SECTION 3. Module level: SFARI genes are not enriched in modules from gene co-expression networks that are 
# strongly correlated with ASD diagnosis

# Fig 5. SFARI gene enrichment in modules does not correlate with ASD diagnosis status
# Size: 4 x 5 in for main figure and 3 x 7 for supplementary figures
modules_dataset = modules_dataset %>% left_join(genes_info %>% dplyr::select(ID, `gene-score`), by = 'ID') %>% 
  filter(Module != 'gray')
# We need the entrez ID of the genes for this
getinfo = c('ensembl_gene_id','entrezgene')
mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='feb2014.archive.ensembl.org')
biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), values=modules_dataset$ID, mart=mart) %>%
                 left_join(modules_dataset %>% dplyr::select(ID,Module,`gene-score`), by = c('ensembl_gene_id'='ID'))
# We need to build a term2gene dataframe with the genes and their SFARI Scores
term2gene = biomart_output %>% mutate('term' = ifelse(`gene-score` == 'Others', 'Others', 'SFARI'), 
                                      'gene' = entrezgene) %>%  dplyr::select(term, gene) %>% distinct
modules = modules_dataset %>% pull(Module) %>% unique %>% as.character
enrichment_data = data.frame('Module' = modules, 'size' = 0, 'pval_ORA' = 0, 'padj_ORA' = 0)
for(i in 1:length(modules)){
  module = modules[i]
  genes_in_module = biomart_output$entrezgene[biomart_output$Module==module]
  ORA_module = enricher(gene = genes_in_module, universe = biomart_output$entrezgene %>% as.character, 
                        pAdjustMethod = 'BH', TERM2GENE = term2gene, 
                        pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize = 50000) %>% 
    data.frame %>% dplyr::select(-geneID,-Description)
  ORA_pval = ifelse('SFARI' %in% ORA_module$ID, ORA_module$pvalue[ORA_module$ID=='SFARI'], 1)
  ORA_padj = ifelse('SFARI' %in% ORA_module$ID, ORA_module$p.adjust[ORA_module$ID=='SFARI'], 1)
  enrichment_data[i,-1] = c(length(genes_in_module), ORA_pval, ORA_padj)
}
enrichment_data = classification_dataset %>% mutate(ID = rownames(.)) %>% dplyr::select(ID, MTcor) %>%
                  inner_join(modules_dataset, by = 'ID') %>% dplyr::select(Module, MTcor) %>% 
                  distinct(MTcor, Module) %>% left_join(enrichment_data, by = 'Module') %>% 
                  mutate('enriched_in_SFARI' = padj_ORA<0.01)

module_diagnosis_vs_SFARI = enrichment_data %>% ggplot(aes(MTcor, 1-pval_ORA, size=size)) + 
                            geom_point(aes(color = enriched_in_SFARI), alpha = 0.5) + 
                            geom_smooth(color = '#cccccc', size = 0.5,alpha = 0.1) + 
                            xlab('Module-Diagnosis correlation') + ylab('Enrichment in SFARI Genes') + 
                            guides(size = FALSE) + labs(color = 'Enrichment is statistically significant') + 
                            theme_minimal() + theme(legend.position = 'bottom')

# Fig 6. SFARI gene enrichment in modules is related to the mean level of expression of the genes in the module
# Size: 4 x 5 in for main figure and 3 x 7 for supplementary figures
mean_expr_by_module = data.frame('ID' = rownames(datExpr), 'meanExpr' =rowMeans(datExpr)) %>%
                      inner_join(modules_dataset %>% dplyr::select(-`gene-score`), by = 'ID') %>% 
                      dplyr::select(-ID) %>% group_by(Module) %>% summarise(meanExpr = mean(meanExpr)) %>% 
                      left_join(enrichment_data, by = 'Module')

mean_expr_vs_SFARI = mean_expr_by_module %>% ggplot(aes(meanExpr, 1-pval_ORA, size=size)) + 
                      geom_point(aes(color = enriched_in_SFARI), alpha = 0.5) + 
                      geom_smooth(color = '#cccccc', size = 0.5,alpha = 0.1) + 
                      xlab('Mean expression') + ylab('Enrichment in SFARI Genes') + 
                      guides(size = FALSE) + labs(color = 'Enrichment is statistically significant') + 
                      theme_minimal() + theme(legend.position = 'bottom')


rm(i, module, genes_in_module, ORA_module, ORA_pval, ORA_padj, getinfo, mart, term2gene)


################################################################################################################
# SECTION 4. Systems level: Whole co-expression network integration facilitates discovery of novel candidate SFARI Genes

# Fig 7. Scatter plot of the relation between the level of expression of the genes and the probabilities assigned 
# by the classification model
# Size: 4 x 5 in for main figure and 3 x 7 for supplementary figures
mean_expr_vs_prob = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>%
                    inner_join(biased_classification_model$predictions, by = 'ID') %>%
                    ggplot(aes(meanExpr, prob)) + geom_point(alpha=0.15, color='#0099cc') +
                    geom_smooth(method='gam', color='gray', alpha=0.2) +
                    xlab('Mean expression') + ylab('Original model probability') +
                    theme_minimal()

# Fig 8. Scatter plot of the genes in the unbiased classification dataset
# Size: 4 x 5 in for main figure and 3 x 7 for supplementary figures
mean_expr_vs_corrected_prob = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>%
                              inner_join(unbiased_classification_model$predictions, by = 'ID') %>%
                              ggplot(aes(meanExpr, prob)) + geom_point(alpha=0.15, color='#0099cc') +
                              geom_smooth(method='gam', color='gray', alpha=0.2) +
                              xlab('Mean expression') + ylab('Unbiased model probability') +
                              theme_minimal()


# Table 1. Performance metrics of the two classification models used as well as a third model using a shuffling 
# of the SFARI labels in the data
pm_table = data.frame('Performance Metrics' = colnames(pm_unbiased_model[,-1]),
                      'Biased Model' = pm_biased_model[,-1] %>% apply(2,mean) %>% round(2) %>% unname,
                      'SD Biased Model' = pm_biased_model[,-1] %>% apply(2,sd) %>% round(4) %>% unname,
                      'Unbiased Model' = pm_unbiased_model[,-1] %>% apply(2,mean) %>% round(2) %>% unname,
                      'SD Unbiased Model' = pm_unbiased_model[,-1] %>% apply(2,sd) %>% round(2) %>% unname,
                      'Random Model' = pm_shuffled_model[,-1] %>% apply(2,mean) %>% round(2) %>% unname,
                      'SD Random Model' = pm_shuffled_model[,-1] %>% apply(2,sd) %>% round(2) %>% unname) %>% 
            transpose %>% dplyr::slice(-1)
colnames(pm_table) = colnames(pm_biased_model[,-1])
rownames(pm_table) = c('BiasedModel', 'SDBiasedModel', 'UnbiasedModel', 'SDUnbiasedModel', 'RandomModel', 'SDRandomModel')

pm_table %>% dplyr::select(AUC, MLP, Balanced_Accuracy) %>% xtable(type = 'latex')

# Table 2. Top 10 non-SFARI genes with the highest probabilities assigned by the unbiased model
top_scored_genes = unbiased_classification_model$predictions %>% filter(SFARI == FALSE) %>% 
                   arrange(desc(prob)) %>% top_n(10, wt = prob) %>% dplyr::rename('Gene' = hgnc_symbol) %>% 
                   dplyr::select(Gene, prob) %>% mutate('Literature Review' = '') 

top_scored_genes %>% xtable(type = 'latex')


################################################################################################################
# SECTION 5. Other ASD scoring systems

all_scores = unbiased_classification_model$predictions %>% dplyr::select(-SFARI) %>% 
             dplyr::rename('gene'=hgnc_symbol, 'our score'=prob) %>%
             full_join(genes_info %>% dplyr::rename('gene' = hgnc_symbol) %>%
                       mutate(SFARI = factor(ifelse(`gene-score` == 'Others', NA, `gene-score`),
                                             levels = c('3','2','1'))), by = 'gene') %>%
             full_join(DisGeNET[['asd']] %>% dplyr::rename('gene' = gene_symbol, 'DisGeNET' = score), by = 'gene') %>% 
             full_join(krishnan, by = 'gene') %>% 
             full_join(TADA, by = 'gene') %>% 
             right_join(data.frame('gene' = genes_info$hgnc_symbol, 'meanExpr' = rowMeans(datExpr)), by = 'gene') %>%
             dplyr::select(gene, `our score`, DisGeNET, SFARI, Krishnan, TADA, meanExpr)


# Fig 9. Pairwise correlation between the different ASD scoring systems
# Size 5 x 5 in and crop
corrplot.mixed(hetcor(all_scores[-c(1,2,7)], use = 'pairwise.complete.obs')$correlations,
               p.mat = cor.mtest(all_scores[,-c(1,2,7)] %>% mutate(SFARI = SFARI %>% as.numeric))$p, sig.level = 0.05,
               lower = 'number', lower.col = '#666666', number.cex = 1, tl.col = '#666666')


# Table 3. Correlation between different ASD scoring systems and the mean level of expression of the genes
cor_score_mean_expr = data.frame('SFARI' = c(hetcor(all_scores$SFARI, all_scores$meanExpr)$correlations[1,-1],
                                             cor.test(all_scores$SFARI %>% as.numeric, all_scores$meanExpr)$p.value),
                                 'Krishnan' = c(hetcor(all_scores$Krishnan, all_scores$meanExpr)$correlations[1,-1],
                                                cor.test(all_scores$Krishnan, all_scores$meanExpr)$p.value),
                                 'TADA' = c(hetcor(all_scores$TADA, all_scores$meanExpr)$correlations[1,-1],
                                            cor.test(all_scores$TADA, all_scores$meanExpr)$p.value),
                                 'DisGeNET' = c(hetcor(all_scores$DisGeNET, all_scores$meanExpr)$correlations[1,-1],
                                                cor.test(all_scores$DisGeNET, all_scores$meanExpr)$p.value))
rownames(cor_score_mean_expr) = c('Correlation','p-value')
cor_score_mean_expr %>% xtable(type = 'latex')

# Table 4. Number of genes associated to different neuronal disorders according to DisGeNET and percentage of 
# genes that belong to the SFARI Genes list
disgenet_info = data.frame('gene_symbol' = genes_info$hgnc_symbol, 'meanExpr' = rowMeans(datExpr),
                           'SFARI' = !genes_info$`gene-score`=='Others') %>% 
                left_join(DisGeNET[['asd']] %>% dplyr::rename('ASD'=score), by = 'gene_symbol') %>%
                left_join(DisGeNET[['scz']] %>% dplyr::rename('Scz'=score), by = 'gene_symbol') %>%
                left_join(DisGeNET[['bd']] %>% dplyr::rename('BD'=score), by = 'gene_symbol') %>%
                left_join(DisGeNET[['id']] %>% dplyr::rename('ID'=score),by='gene_symbol')%>%
                left_join(DisGeNET[['dd']] %>% dplyr::rename('DD'=score), by='gene_symbol') %>%
                left_join(DisGeNET[['ai']] %>% dplyr::rename('CAI'=score), by = 'gene_symbol') %>%
                dplyr::select(gene_symbol, meanExpr, ASD, Scz, BD, ID, DD, CAI, SFARI)

perc_SFARI = data.frame('SFARI' = disgenet_info %>% filter(SFARI == TRUE) %>% 
                                  dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                  apply(2, function(x) sum(!is.na(x))),
                        'Total' = disgenet_info %>% dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                  apply(2, function(x) sum(!is.na(x)))) %>%
             mutate(SFARI = 100*SFARI/Total)
rownames(perc_SFARI) = c('ASD','Scz','BD','ID','DD','CAI')

perc_SFARI %>% t %>% xtable(type = 'latex')

# Fig 10. Box plots of gene scores from different disorders comparing SFARI Genes with the rest of the genes
# Size: 3 x 6.5 in
SFARI_in_other_disorders = disgenet_info %>% dplyr::select(-gene_symbol) %>% 
                           reshape2::melt(id.vars = c('meanExpr', 'SFARI')) %>% 
                           dplyr::rename('Disorder' = variable) %>% filter(!is.na(value)) %>%
                           ggplot(aes(SFARI, value, fill = SFARI)) + 
                           geom_boxplot(outlier.colour='gray', outlier.shape='o', outlier.size=3) + 
                           facet_grid(~Disorder) + scale_y_log10(limits = c(NA, 0.7)) +
                           stat_compare_means(label = 'p.signif', method = 't.test', method.args = list(var.equal=F)) +
                           xlab('') + ylab('DisGeNET Score') + theme_minimal() + theme(legend.position = 'bottom')

# Table 5. Correlation between the scores associated to different disorders by DisGeNET and the mean level of 
# expression of the genes
cor_score_mean_expr = data.frame('Genes' =  disgenet_info %>% dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                           apply(2, function(x) sum(!is.na(x))), 
                                 'Correlation' = disgenet_info %>% dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                                 apply(2, function(x) cor.test(x, disgenet_info$meanExpr)$estimate),
                                 'p-value' = disgenet_info %>% dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                             apply(2, function(x) cor.test(x, disgenet_info$meanExpr)$p.value))

# Table 6. Correlation between the scores associated to different disorders by DisGeNET, removing the SFARI genes, 
# and the mean level of expression of the genes
cor_score_mean_expr = data.frame('Genes' =  disgenet_info %>% filter(SFARI == FALSE) %>% 
                                            dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                            apply(2, function(x) sum(!is.na(x))),
                                 'Correlation' = disgenet_info %>% filter(SFARI == FALSE) %>%
                                                 dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                                 apply(2, function(x) 
                                                   cor.test(x, disgenet_info$meanExpr[!disgenet_info$SFARI])$estimate),
                                 'p-value' = disgenet_info %>% filter(SFARI == FALSE) %>%
                                             dplyr::select('ASD','Scz','BD','ID','DD','CAI') %>% 
                                             apply(2, function(x) 
                                               cor.test(x, disgenet_info$meanExpr[!disgenet_info$SFARI])$p.value))

cor_score_mean_expr %>% t %>% xtable(type = 'latex')


################################################################################################################
################################################################################################################
# PLOTS FOR THE METHODS SECTION

# Fig 11. PCA plots of samples before and after preprocessing
# Size: 4 x 5 in for main figure and 3 x 7 for supplementary figures
pca = datExpr %>% t %>% prcomp
pca_samples = data.frame( 'PC1' = pca$x[,1], 'PC2' = pca$x[,2], 'Diagnosis'=datMeta$Diagnosis) %>% 
  ggplot(aes(PC1, PC2, color = Diagnosis, shape = Diagnosis)) + 
  geom_point(alpha=0.9) + coord_fixed() + theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste0('PC1 (',round(100*summary(pca)$importance[2,1],1),'%)')) +
  ylab(paste0('PC2 (',round(100*summary(pca)$importance[2,2],1),'%)'))

#ggarrange(pca_samples_Gupta, pca_samples_Wright, nrow=1, common.legend = TRUE, legend='bottom', widths = c(0.59, 0.41))

rm(pca)

# Fig 12. Bias correction
# Size: 3.5 x 7 in for main figure and 7 x 7 for supplementary figures

# Fig 12A. Balanced accuracy and bias through each of the iterations of the bias correction algorithm
plot_info = data.frame('iter' = 1:length(unbiased_classification_model$BiasVector), 
                       'Bias' = unbiased_classification_model$BiasVector, 
                       'BalancedAccuracy' = unbiased_classification_model$BalancedAccuracyVector) %>% 
  reshape2::melt(id.vars = 'iter') %>% 
  mutate(variable = ifelse(variable == 'BalancedAccuracy','Balanced accuracy', 'Bias'))

bias_correction_iters = plot_info %>% ggplot(aes(x=iter, y=value, color = variable)) + geom_line() + 
  xlab('Iteration') + ylab('Value') + theme_minimal() + labs(color = 'Metric ') +
  theme(legend.position = 'bottom')


# Fig 12B. Weights assigned to each gene by the final iteration of the bias correction algorithm based on their 
# mean expression and label
mean_expr = data.frame('ID' = rownames(datExpr), 'meanExpr' = rowMeans(datExpr)) %>%
  left_join(unbiased_classification_model$predictions, by = 'ID') %>% filter(n>0) %>%
  mutate('meanExpr_std' = (meanExpr-mean(meanExpr))/sd(meanExpr))

w_hat = exp(unbiased_classification_model$lambda*mean_expr$meanExpr_std) # inverso a mean expr
w0 = 1/(1+w_hat) # prop a mean expr
w = 1/(1+w_hat)

# inv mean expr Positives, prop Negatives
w[mean_expr$SFARI %>% as.logical] = w[mean_expr$SFARI %>% as.logical]*w_hat[mean_expr$SFARI %>% as.logical]
plot_data = data.frame('meanExpr' = mean_expr$meanExpr, 'w_hat' = w_hat, 'w0' = w0, 'w' = w, 
                       'SFARI' = mean_expr$SFARI, 'pred' = mean_expr$pred)

bias_correction_weights = plot_data %>% ggplot(aes(meanExpr, w, color = SFARI)) + geom_point(alpha = 0.2) + 
  xlab('Mean expression') + ylab('Weight') + theme_minimal() + theme(legend.position='bottom')


#grid.arrange(bias_correction_iters, bias_correction_weights, nrow = 1, widths = c(0.45, 0.55))
