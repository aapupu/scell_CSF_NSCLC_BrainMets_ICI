import loompy as lp
# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "Mono_scenic.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "Mono.h5ad"

# path to pyscenic output
f_pyscenic_output = "Mono_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = 'Mono_scenic_integrated-output.loom'

row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( adata.obs['nFeature_RNA']) ,
    "nUMI": np.array(  adata.obs['nCount_RNA']) ,
}

lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs )
f_tfs = "pyscenic_file/allTFs_hg38.txt" 
!pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 10

adjacencies = pd.read_csv("adj.csv", index_col=False, sep='\,')

import glob
# ranking databases
f_db_glob = "pyscenic_file/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "pyscenic_file/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

!pyscenic ctx adj.csv \
    {f_db_names} \
    --annotations_fname {f_motif_path} \
    --expression_mtx_fname {f_loom_path_scenic} \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 10

!pyscenic aucell \
    {f_loom_path_scenic} \
    reg.csv \
    --output {f_pyscenic_output} \
    --num_workers 10

lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

from pyscenic.rss import regulon_specificity_scores
rss_cellType = regulon_specificity_scores( auc_mtx, adata.obs['Celltype'] )
plot_rss(rss_cellType, c, top_n=10, max_n=None)

