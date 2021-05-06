import numpy as np
import scanpy as sc
import pandas as pd
import os
import pandas as pd
import gc
import dask.dataframe

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300,dpi_save=600, facecolor='white')

adata_tmp=sc.read_csv("sct_dge.csv")
dge=pd.DataFrame(adata_tmp.X)
dge.index=pd.read_csv("sct_gene.csv")["x"]
dge.columns=pd.read_csv("sct_cell.csv")["x"]
adata=sc.AnnData(dge.T)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.tsne(adata, n_pcs=50)
sc.tl.leiden(adata, resolution=2.5)#2~3
sc.pl.tsne(adata, color=['leiden'],legend_fontsize=5,save="leiden.pdf")

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges','scores', 'pvals_adj']}).to_csv("marker.csv")

adata1.write('sct_atlas.h5ad',compression=True)