import scanpy as sc
import pandas as pd
import os
results = pd.read_csv("/lvdata/wzb/pipline/function/Cellphone/mouse2human.xls",index_col=0,sep="\t")
adata = sc.read_h5ad("/lvdata/wzb/scRNA/FW2024-021/GEX/result.h5ad")
out_dir = "/lvdata/wzb/scRNA/FW2024-021/2024.05.30/"
# print(results)
def trans(adata,results,out_dir = os.getcwd()):
    adata1 = adata.raw.to_adata()
    adata1 = adata1[:,adata1.var_names.isin(results.index.to_list())]
    tmp = adata1.var.merge(results,left_index=True, right_index=True,how='left',validate ="1:1")
    # print(tmp)
    tmp = tmp[~tmp['human_gene'].duplicated(keep='first')]

    adata1 = adata1[:,adata1.var_names.isin(tmp.index.to_list())]
    adata1.var = tmp
    adata1.var = adata1.var.set_index('human_gene')
    # adata1.obs['cell_type'] = adata1.obs['cell_type_1']
    # sc.pl.umap(adata1,color='cell_type',show=False)
    adata1.raw = adata1
    adata1.write_h5ad(os.path.join(out_dir,'human.h5ad'))

trans(adata,results,out_dir)