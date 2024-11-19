RelRankEnrich <- function(exp, genelist){
  if (!is.matrix(exp) || is.null(rownames(exp)) || length(genelist) == 0) {
    stop("Invalid input. 'exp' must be a non-empty matrix and 'genelist' must be a non-empty list.")
  }
  row_names = rownames(exp)
  num_genes = nrow(exp)
  num_samples = ncol(exp)
  
  R = matrixStats::rowRanks(exp, preserveShape = T, ties.method = 'average')
  
  result_matrix = matrix(0, nrow = length(genelist), ncol = num_samples)
  rownames(result_matrix) <- names(genelist)
  colnames(result_matrix) <- colnames(exp)
  
  for (i in names(genelist)) {
    gene_set = genelist[[i]]
    valid_genes = intersect(row_names, gene_set)
    
    if (length(valid_genes) == 0) {
      warning(paste("No valid genes for gene set", i, ". Setting enrichment score to 0."))
      es <-  rep(0,num_samples)
      
    }else if(length(valid_genes) == 1){
      gene_set_idx <- match(valid_genes,row_names)
      es <- R[gene_set_idx,]
      es <- es / num_samples
      
    }else{
      gene_set_idx <- match(valid_genes,row_names)
      es <- colSums(R[gene_set_idx,])
      es <- es / length(valid_genes)
      es <- es / num_samples
    }
    result_matrix[i, ] = es
  }
  return(result_matrix)
}

CD4_Trm_CXCR6 <- c('CD3D','CD2','CD5','CD4',
                   'CD6','CD28','HMGB2','SLAMF6',
                   'NFATC2','NFATC3',
                   'CD44','WNK1','CXCR6','RGS1',
                   'PD-1','TIM-3',
                   'EOMES','PRDM1',
                   'IFN-gamma','BLC','IL-21',
                   'GZMH','CST7',
                   'SKAP1','HSP40','P2RY10','SH3KBP1')
genelist <- list(CD4_Trm_CXCR6=CD4_Trm_CXCR6)
es <- RelRankEnrich(protein_exp, genelist) %>% as.data.frame() %>% t()
protein_meta$`RelRankEnrich of CD4_Trm_CXCR6` <- es[,1]
protein_meta$Group <- factor(protein_meta$Group,levels = c("CR_pre",'PR_pre','NR_pre','NR_post'))
ggboxplot(data =protein_meta,x = "Group",y = "RelRankEnrich of CD4_Trm_CXCR6"  , color ="Group",
          add = "jitter",palette = c("#00AFBB",'orange', "#FC4E07","#3C5488FF" ))+
  xlab("")+ylab('RelRankEnrich of protein list')+ggtitle('CD4_Trm_CXCR6')+
  theme(plot.title = element_text(hjust = 0.5,size = 13))+
  stat_compare_means(comparisons = list(c("CR_pre",'NR_pre'),c('NR_pre','NR_post'),c('PR_pre','NR_pre'),c("CR_pre",'PR_pre')),
                     method = 't.test',size = 4)+NoLegend()


ProteinType=c(rep('Marker of T cell',4),rep('Regulation and activation',6),rep('Migration and residence',4), 
           rep('Immune Checkpoint',2),rep('Transcription factor',2),rep('Cytokine',3),rep('Cytotoxicity',2),
           rep('Others',4))
names(ProteinType)=CD4_Trm_CXCR6
ProteinType=as.data.frame(ProteinType)
ProteinType$ProteinType <- factor(ProteinType$ProteinType,levels = c('Marker of T cell','Regulation and activation','Migration and residence',
                                                                     'Immune Checkpoint', 'Transcription factor','Cytokine','Cytotoxicity','Others'))
ann_colors = list(
  ProteinType = c(`Marker of T cell` = "#66C2A5",`Regulation and activation` = "#8DA0CB",`Migration and residence` = "#E78AC3",
                  `Immune Checkpoint` = "#A6D854",`Transcription factor` = "#FFD92F",`Cytokine` = "#E5C494",
                  `Cytotoxicity` = "#B3B3B3",Others=  "#D9D9D9"))
protein_exp_mean_sub <- protein_exp_mean[,CD4_Trm_CXCR6]
pheatmap(t(protein_exp_mean_sub),scale = 'row', color = colours,annotation_row = ProteinType, annotation_colors = ann_colors,
        cluster_cols = F,cluster_rows = F,border_color = NA,angle_col = '45'
        )



