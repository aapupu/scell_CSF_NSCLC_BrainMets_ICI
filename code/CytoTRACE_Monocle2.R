# CytoTRACE
cytotraceresult <- CytoTRACE::CytoTRACE(mat = as.matrix(cDC@assays$RNA@counts),ncores = 20,
                                        enableFast = F, batch=as.character(cDC$Batch))

# Monocle2
diff_test_res <- differentialGeneTest(monocle_data,
                                      fullModelFormulaStr = "~Major_Celltype",cores = 20)
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
ordering_genes <- rownames(diff_test_res[1:300,])
monocle_data_select <- setOrderingFilter(monocle_data, ordering_genes)
monocle_data_select <- reduceDimension(monocle_data_select, max_components = 2,
                                       reduction_method="DDRTree",residualModelFormulaStr="~Batch")
monocle_data_select <- orderCells(monocle_data_select)


diff_test_CR <- differentialGeneTest(monocle_data_select[,monocle_data_select$Major_Celltype %in% c('cDC2','migDC') & 
                                                             monocle_data_select$Group5 =='CR_pre'],
                                       fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 30)
diff_test_CR <- diff_test_CR[order(diff_test_CR$qval),]
sig_gene_names1 <- diff_test_CR$gene_short_name[diff_test_CR$qval<0.01]
sig_gene_names1
