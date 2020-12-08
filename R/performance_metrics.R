
setwd('~/PhD/Paper/')

library(doParallel)

# Load necessary scripts
source('R/packages.R')
source('R/auxiliary_functions.R')
source('R/main_functions.R')

#loadd(Gandal_dataset, classification_dataset, new_SFARI_dataset)
load('Results/Gandal_dataset.RData')
load('Results/Gandal_dataset.RData')
load('Results/Gandal_dataset.RData')

n_iter = 100
seeds = 1:n_iter*1000 # We generate 100 seeds inside the model from the original one, this way they don't overlap between iterations

# Setup parallel backend to use many processors
cores = detectCores()
cl = makeCluster(cores-2)
registerDoParallel(cl)

# Create matrix with the performance metrics from each iteration
performance_metrics = foreach(i=seeds, .combine=rbind, .errorhandling = 'remove') %dopar% {
  
  source('R/packages.R')
  library(dplyr)
  
  ###################################################################################
  # Shuffle SFARI labels
  # set.seed(i)
  # random_classification_dataset = classification_dataset
  # random_classification_dataset$SFARI = sample(classification_dataset$SFARI)
  # random_SFARI_dataset = new_SFARI_dataset[new_SFARI_dataset$ID %in% rownames(random_classification_dataset),]
  # random_SFARI_dataset$ID[!is.na(random_SFARI_dataset$`gene-score`)] =
  #   rownames(random_classification_dataset[random_classification_dataset$SFARI,])
  # random_SFARI_dataset = random_SFARI_dataset %>% distinct(ID, .keep_all = TRUE)
  # ###################################################################################

  results = classification_model(Gandal_dataset, classification_dataset, SFARI_dataset, i, FALSE)

  temp = as.vector(c('seed' = i, unlist(results$pm_final_model)))
  
}

# Stop clusters
stopCluster(cl)

colnames(performance_metrics) = c('Seed','Accuracy','Precision','Recall','F1','AUC','MLP','Balanced_Accuracy')

write.csv(performance_metrics, file = 'Results/performance_metrics_biased_model.csv', row.names=FALSE)
#write.csv(performance_metrics, file = 'Results/performance_metrics_shuffled_labels_model_8.csv', row.names=FALSE)


################################################################################################################
################################################################################################################
# # PLOT RESULTS
# 
# # performance_metrics = read.csv('Results/performance_metrics_final_model.csv')
# 
# # colMeans(performance_metrics)
# # apply(performance_metrics, 2, sd)
# 
# # top_genes = nbiased_classification_model$predictions %>% filter(!SFARI) %>% 
# # dplyr::arrange(desc(prob)) %>% dplyr::select(ID, hgnc_symbol, prob) %>% 
# # left_join(data.frame(ID = rownames(classification_dataset), MTcor = classification_dataset$MTcor, 
# # GS = classification_dataset$GS), by = 'ID') %>% dplyr::select(hgnc_symbol, GS, MTcor, prob) %>% 
# # mutate(LitReview = '') %>% head(10)
# # xtable(top_genes, type = 'latex')
# 
# # MEAN EXPRESSION
# # PCA
# library(viridis)
# library(ggpubr)
# library(rstatix)
# library(gridExtra)
# 
# SFARI_colour_hue = function(r) {
#   pal = c('#FF7631','#FFB100','#E8E328','#8CC83F','#62CCA6','#59B9C9','#b3b3b3','#808080','gray','#d9d9d9')[r]
# }
# 
# pca = datExpr %>% prcomp
# plot_data = data.frame( 'PC1' = pca$x[,1], 'PC2' = pca$x[,2], 'MeanExpr'=rowMeans(datExpr))
# plot_data %>% ggplot(aes(PC1, PC2, color=MeanExpr)) + geom_point(alpha=0.3) + theme_minimal() + 
#               scale_color_viridis() + ggtitle('PCA of Genes') +
#               xlab(paste0('PC1 (',round(100*summary(pca)$importance[2,1],1),'%)')) +
#               ylab(paste0('PC2 (',round(100*summary(pca)$importance[2,2],1),'%)'))
# 
# # over vs under-expressed
# plot_data = data.frame(ID = rownames(datExpr), MeanExpr = rowMeans(datExpr)) %>% 
#   left_join(genes_info, by = 'ID') %>%
#   mutate(Group = ifelse(log2FoldChange>0, 'LFC>0', 'LFC<0')) %>%
#   mutate(Group = factor(Group, levels = c('LFC<0', 'LFC>0')))
# p1 = plot_data %>% ggplot(aes(Group, MeanExpr, fill=Group)) + ylab('Mean Expression') + xlab('Group') +
#   geom_boxplot(outlier.colour='gray', outlier.shape='o', outlier.size=3) + 
#   stat_compare_means(label = 'p.signif', method = 't.test', method.args = list(var.equal = FALSE)) +
#   theme_minimal() + theme(legend.position='none')
# 
# wt = plot_data %>% wilcox_test(MeanExpr~gene.score, p.adjust.method='BH') %>% add_x_position(x = 'group')
# increase = 1
# base = 15.5
# pos_y_comparisons = c(base + c(0,1,3,5)*increase, base+c(0,2,4)*increase, base+0:1*increase, base)
# p2 = plot_data %>% ggplot(aes(gene.score, MeanExpr)) + 
#   geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, aes(fill=gene.score)) +
#   stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .01) +
#   scale_fill_manual(values=SFARI_colour_hue(r=c(1:3,8,7))) + 
#   xlab('SFARI Gene Scores') + ylab('Mean Expression') + 
#   theme_minimal() + theme(legend.position='none')
# 
# grid.arrange(p1, p2, nrow=1, widths = c(0.36, 0.65))
# 
# # Bias in model
# loadd(biased_classification_model)
# p1 = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>% 
#      left_join(biased_classification_model$predictions, by = 'ID') %>% 
#      ggplot(aes(meanExpr, prob)) + geom_point(alpha=0.1, color='#0099cc') +
#      geom_smooth(method='gam', color='gray', alpha=0.2) + 
#      xlab('Mean Expression') + ylab('Original Probability') +
#      theme_minimal() + ggtitle('Original model')
# 
# p2 = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>% 
#      left_join(unbiased_classification_model$predictions, by = 'ID') %>% 
#      ggplot(aes(meanExpr, prob)) + geom_point(alpha=0.1, color='#ff9900') +
#      geom_smooth(method='gam', color='gray', alpha=0.2) + 
#      xlab('Mean Expression') + ylab('Corrected Probability') +
#      theme_minimal() + ggtitle('Bias corrected model')
# 
# grid.arrange(p1, p2, nrow=1)
# 
# # LFC
# 
# wt = genes_info %>% mutate(abs_lfc = abs(log2FoldChange)) %>%
#   wilcox_test(abs_lfc~gene.score, p.adjust.method = 'BH') %>% add_x_position(x = 'group')
# 
# increase = 0.05
# base = 0.45
# pos_y_comparisons = c(base + c(0,1,3,5)*increase, base+c(0,2,4)*increase, base+0:1*increase, base)
# genes_info %>% ggplot(aes(gene.score, abs(log2FoldChange))) + 
#   geom_boxplot(outlier.colour='#cccccc', outlier.shape='o', outlier.size=3, aes(fill=gene.score)) +
#   stat_pvalue_manual(wt, y.position = pos_y_comparisons, tip.length = .005) +
#   scale_fill_manual(values=SFARI_colour_hue(r=c(1:3,8,7))) + 
#   coord_cartesian(ylim = c(0, max(pos_y_comparisons))) +
#   xlab('SFARI Gene Scores') + ylab('LFC Magnitude') + 
#   theme_minimal() + theme(legend.position='none')
# 
# p1 = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>% 
#   left_join(unbiased_classification_model$predictions, by = 'ID') %>% 
#   left_join(genes_info, by = 'ID') %>% 
#   mutate(Group = ifelse(log2FoldChange>0, 'LFC>0', 'LFC<0')) %>%
#   mutate(Group = factor(Group, levels = c('LFC<0', 'LFC>0'))) %>%
#   filter(Group == 'LFC<0') %>%
#   ggplot(aes(abs(log2FoldChange), prob)) + geom_point(alpha=0.15, color = '#F9766E') +
#   geom_smooth(method='lm', color='gray', alpha=0.2) + scale_x_sqrt() +
#   xlab('LFC Magnitude') + ylab('Corrected Probability') +
#   theme_minimal() + theme(legend.position = 'none') + ggtitle('LFC<0 Genes')
# 
# p2 = data.frame(ID = rownames(datExpr), meanExpr = rowMeans(datExpr)) %>% 
#   left_join(unbiased_classification_model$predictions, by = 'ID') %>% 
#   left_join(genes_info, by = 'ID') %>% 
#   mutate(Group = ifelse(log2FoldChange>0, 'LFC>0', 'LFC<0')) %>%
#   mutate(Group = factor(Group, levels = c('LFC<0', 'LFC>0'))) %>%
#   filter(Group == 'LFC>0') %>%
#   ggplot(aes(log2FoldChange, prob)) + geom_point(alpha=0.15, color = '#00BFC4') +
#   geom_smooth(method='lm', color='gray', alpha=0.2) + scale_x_sqrt() +
#   xlab('LFC Magnitude') + ylab('Corrected Probability') +
#   theme_minimal() + theme(legend.position = 'none') + ggtitle('LFC>0 Genes')
# 
# grid.arrange(p1, p2, nrow=1)

