import scanpy as sc
import scvelo as scv
# CD4_T
sc.tl.paga(adata, groups='Celltype')
sc.tl.paga(adata_CR_pre, groups='Celltype')
sc.tl.paga(adata_NR_pre, groups='Celltype')
sc.tl.paga(adata_NR_post, groups='Celltype')

adata.uns['iroot'] = np.flatnonzero(adata.obs['Celltype']  == 'CD4_Tn_CCR7')[0]
sc.tl.dpt(adata)
sc.tl.draw_graph(adata,init_pos='paga')

# Monocyte & Macropahge
ldata = scv.read('merge.loom', cache=True)
barcodes = [bc.replace(':','_') for bc in ldata.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] for bc in barcodes]
ldata.obs.index = barcodes
ldata = ldata[adata.obs.index, :]
adata = scv.utils.merge(adata, ldata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=15, n_neighbors=30)
scv.tl.recover_dynamics(adata,n_jobs=20)
scv.tl.velocity(adata,mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=1, node_size_scale=1,title="",legend_loc='right',
            palette=["#725663FF","#C71000FF", "#008EA0FF","#8A4198FF" ,"#5A9599FF", "#FF6348FF" ,"#84D7E1FF" ,"#FF95A8FF"],
            save='Mac_Mono_velocity_paga.png')
