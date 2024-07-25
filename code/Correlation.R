# CSF
df1 <- data.frame(prop.table(table(sce$Celltype,sce$Sample),margin = 2))
colnames(df1) <- c("Celltype","Sample","Proportion")
df1_wide<-dcast(df1,Sample~df1$Celltype,value.var = 'Proportion')
df1_corr <- cor(df1_wide[,-1],method = 'pearson') %>% as.data.frame()
ComplexHeatmap::pheatmap(df1_corr,border_color = NA,fontsize = 7,cellwidth = 6,cellheight = 6,
                   breaks = seq(-1,1,length.out = 101), display_numbers = F,show_colnames=F,
                   treeheight_row=F,treeheight_col=40,
                   cutree_cols=5,cutree_rows=5,border=NA,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean", clustering_method = "complete",
                   color =colorRampPalette(rev(brewer.pal(n =11, name ="RdBu")))(100),name='PCC',
                   top_annotation = HeatmapAnnotation(
                     foo = anno_block(
                       gp = gpar(fill = c(4,3,6,7,2)),
                       labels = c("Cluster1", "Cluster2", "Cluster3","Cluster4","Cluster5"), 
                       labels_gp = gpar(col = "black", fontsize = 8))
                   ), right_annotation = rowAnnotation(
                     foo = anno_block(
                       gp = gpar(fill = c(4,3,6,7,2)),
                       #labels = c("Cluster1", "Cluster2", "Cluster3","Cluster4","Cluster5"), 
                       labels_gp = gpar(col = NA))
                   ))

# CSF_tissue
Celltype_Corr <- function(seurat_obj, assays='RNA',Celltype, method){
  av <- AverageExpression(object = seurat_obj,group.by = Celltype,
                          assays = assays,slot = 'data')[[1]]
  cg <-  names(tail(sort(apply(av,1,sd)),2000))
  
  cor_matrix <- cor(av[cg,],method = method)
  pheatmap(cor_matrix,border_color = NA,
           #breaks = seq(-1,1,length.out = 101), 
           color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))
  return(cor_matrix)
}
CD4_merge <- merge(subset(CD4.combined_sub,Celltype != 'CD4_Treg'),subset(Tissue_sce,Celltype == 'CD4'))
CD4_cor_matrix <- Celltype_Corr(CD4_merge,Celltype = 'Celltype',method = 'spearman',assays = 'RNA')
