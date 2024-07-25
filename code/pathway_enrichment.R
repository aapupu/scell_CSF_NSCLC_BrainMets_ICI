# GO BP
GO_result  <-   enrichGO(gene          = CD4_degs$gene[CD4_degs$cluster=='CD4_Trm_CXCR6'],
                         #universe     = row.names(dge.celltype),
                         OrgDb         = 'org.Hs.eg.db',
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.1)
GO_result <- as.data.frame(GO_result)
GO_result$`-log10P` <- -log10(GO_result$p.adjust)
GO_result <- GO_result[order(GO_result$`-log10P`,decreasing = T),]
GO_result$Description <- factor(GO_result$Description,levels = rev(GO_result$Description)) 
ggplot(data = GO_result[1:10,],mapping = aes(x=Description,y=`-log10P`))+
  geom_bar(stat = "identity",width = 0.85,color="darkred",fill="darkred",)+
  coord_flip()+theme_bw()+
  theme(axis.text.y = element_text(size = 14,color = "black"),
        axis.text.x = element_text(size = 14,color = "black"),
        axis.title = element_text(size=14,color = "black"))+
  ylab("-log10(p.adjust)")+xlab("GO BP")

# KEGG

