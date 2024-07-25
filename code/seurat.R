# scRNA-seq data processing and clustering
sce <- CreateSeuratObject(counts = counts, project = 'ICI', min.cells = 3, min.features = 200) 
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce <- subset(sce, subset = percent.mt < 25 & nFeature_RNA > 500 & nFeature_RNA < 7000)

s_qc_list <- SplitObject(sce,split.by = "Batch")
s_qc_list <- lapply(X = s_qc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list =s_qc_list)

s_qc.anchors <- FindIntegrationAnchors(object.list =s_qc_list, anchor.features = features)
s_qc.combined <- IntegrateData(anchorset = s_qc.anchors)

DefaultAssay(s_qc.combined) <- "integrated"
s_qc.combined <- ScaleData(s_qc.combined, verbose = FALSE, vars.to.regress =c("nFeature_RNA", "nCount_RNA", "percent.mt"))
s_qc.combined <- RunPCA(s_qc.combined, npcs = 30, verbose = FALSE)

s_qc.combined <- RunTSNE(s_qc.combined, reduction = "pca", dims = 1:30)
# s_qc.combined <- RunUMAP(s_qc.combined, reduction = "pca", dims = 1:30)
s_qc.combined <- FindNeighbors(s_qc.combined, reduction = "pca", dims = 1:30)
s_qc.combined <- FindClusters(s_qc.combined)

# AddModuleScore example
sce <- AddModuleScore(sce,features = list(c('CCL2', 'CCL3', 'CCL4','CCL5', 'CXCL10', 'CXCL9', 'IL1B', 'IL6', 'IL7', 'IL15','IL18')),
                     assay = 'RNA',name = "Inflammatory score")
