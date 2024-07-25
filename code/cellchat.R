library(CellChat)

data.input =GetAssayData(object =s_qc.combined_sub,slot = 'data')
cell.use_CR = rownames(meta)[meta$Group5  == "CR_pre"] # extract the cell names from disease data
cell.use_NR = rownames(meta)[meta$Group5  == "NR_pre"] # extract the cell names from disease data
cell.use_NR_post = rownames(meta)[meta$Group5  == "NR_post"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
data.input_CR = data.input[, cell.use_CR]
meta_CR = meta[cell.use_CR, ]        

data.input_NR = data.input[, cell.use_NR]
meta_NR = meta[cell.use_NR, ]   

data.input_NR_post = data.input[, cell.use_NR_post]
meta_NR_post = meta[cell.use_NR_post, ]   

cellchat_CR <- createCellChat(object = data.input_CR, meta = meta_CR, group.by = "Major_Celltype")
cellchat_NR <- createCellChat(object = data.input_NR, meta = meta_NR, group.by = "Major_Celltype")
cellchat_NR_post <- createCellChat(object = data.input_NR_post, meta = meta_NR_post, group.by = "Major_Celltype")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling",'Cell-Cell Contact','ECM-Receptor'))
cellchat_CR@DB <- CellChatDB.use
cellchat_NR@DB <- CellChatDB.use
cellchat_NR_post@DB <- CellChatDB.use

cellchat_CR <- subsetData(cellchat_CR) 
cellchat_NR <- subsetData(cellchat_NR) 
cellchat_NR_post <- subsetData(cellchat_NR_post)

future::plan("multiprocess", workers = 30)

cellchat_CR <- identifyOverExpressedGenes(cellchat_CR)
cellchat_CR <- identifyOverExpressedInteractions(cellchat_CR)

cellchat_NR <- identifyOverExpressedGenes(cellchat_NR)
cellchat_NR <- identifyOverExpressedInteractions(cellchat_NR)

cellchat_NR_post <- identifyOverExpressedGenes(cellchat_NR_post)
cellchat_NR_post <- identifyOverExpressedInteractions(cellchat_NR_post)


cellchat_CR <- computeCommunProb(cellchat_CR)
cellchat_CR <- filterCommunication(cellchat_CR, min.cells = 10)

cellchat_NR <- computeCommunProb(cellchat_NR)
cellchat_NR <- filterCommunication(cellchat_NR, min.cells = 10)

cellchat_NR_post <- computeCommunProb(cellchat_NR_post)
cellchat_NR_post <- filterCommunication(cellchat_NR_post, min.cells = 10)
#
cellchat_CR <- computeCommunProbPathway(cellchat_CR)
cellchat_NR <- computeCommunProbPathway(cellchat_NR)
cellchat_NR_post <- computeCommunProbPathway(cellchat_NR_post)
#
cellchat_CR <- aggregateNet(cellchat_CR)
cellchat_NR <- aggregateNet(cellchat_NR)
cellchat_NR_post <- aggregateNet(cellchat_NR_post)

cellchat_CR <- netAnalysis_computeCentrality(cellchat_CR, slot.name = "netP")
cellchat_NR <- netAnalysis_computeCentrality(cellchat_NR, slot.name = "netP")
cellchat_NR_post <- netAnalysis_computeCentrality(cellchat_NR_post, slot.name = "netP")






