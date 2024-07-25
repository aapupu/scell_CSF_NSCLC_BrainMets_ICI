# GO BP
library(clusterProfiler)
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
library(GSVA)
library(limma)

de_gsva  <- function(exprSet,meta,compare = NULL){
  
  
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}

Group_GSVA <- function(seurat_obj, Group, category='KEGG',compare = compare) {
  expr=as.matrix(seurat_obj@assays$RNA@data)
  
  if(category == 'GO.BP'){                          
    msgdC5 = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
    GO.BPSet = msgdC5 %>% split(x = .$gene_symbol, f = .$gs_name)
    sc_GSVA <- gsva(expr, gset.idx.list = GO.BPSet, kcdf="Gaussian",method = "gsva",
                    parallel.sz=30)   
  }
  if(category == 'KEGG'){
    kegg_df <- read.csv("KEGGREST_WithGene.csv",row.names = 1)
    gene_list <- strsplit(kegg_df$hgnc_symbol,split = ',')
    names(gene_list) <- kegg_df$pathway_name
    sc_GSVA <- gsva(expr, gset.idx.list =gene_list, kcdf="Gaussian",method = "gsva",
                    parallel.sz=30)   
  }
  
  
  #
  meta <- seurat_obj@meta.data[,Group]
  meta <- str_replace(meta,pattern = '-',replacement = '_')
  Diff =de_gsva(exprSet = sc_GSVA,meta = meta,compare = compare)
  Padj_threshold=0.05
  idiff <-Diff[[compare]]
  df <- data.frame(ID = rownames(idiff), score = idiff$t )
  df$group =sapply(1:nrow(idiff),function(x){
    if(idiff[x,"logFC"]>0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("up")}
    else if(idiff[x,"logFC"]<0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("down")}
    else{return("noSig")}
  })
  
  df$hjust = ifelse(df$score>0,1,0)
  df$nudge_y = ifelse(df$score>0,-0.1,0.1)
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  limt = max(abs(df$score))
  Mono_sortdf_sub <- sortdf
  p <- ggplot(Mono_sortdf_sub, aes(ID, score,fill=group)) + 
    geom_bar(stat = 'identity',alpha = 0.7) + 
    scale_fill_manual(breaks=c("down","noSig","up"),
                      values = c("#008020","grey","#08519C"))+
    geom_text(data = Mono_sortdf_sub, aes(label = Mono_sortdf_sub$ID, y = Mono_sortdf_sub$nudge_y),
              nudge_x =0,nudge_y =0,hjust =Mono_sortdf_sub$hjust,
              size = 3)+
    labs(x = paste0('KEGG'," pathways"),
         y=paste0("t value of GSVA score\n",'NR_post-NR_pre'),
         title = 'Monocyte')+
    scale_y_continuous(limits=c(-25,25))+
    coord_flip() + 
    theme_bw() + 
    theme(panel.grid =element_blank())+
    theme(panel.border = element_rect(size = 0.6)
          #panel.border = element_blank()
    )+
    theme(plot.title = element_text(hjust = 0.5,size = 18),
          axis.text.y = element_blank(),
          axis.title = element_text(hjust = 0.5,size = 18),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = limt
    )
  return(list(sortdf,p))
}
B_kegg1 <- Group_GSVA(seurat_obj = B[,B$Group5 %in% c("CR_pre",'NR_pre')],Group = 'Group5',category = "KEGG",compare = 'CR_pre-NR_pre')
