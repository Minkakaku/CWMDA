import os, sys
import loompy as lp
import numpy as np
import anndata as ad
output = "/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/SCENIC"
adata = ad.read_h5ad( "/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/Result.h5ad")
row_attrs = {"Gene": np.array(adata.raw.var_names),}
col_attrs = {"CellID": np.array(adata.raw.obs_names)}
lp.create(os.path.join(output,"sce.loom"),adata.raw.X.transpose(),row_attrs,col_attrs)

os.system(f'sh /lvdata/wzb/pipline/function/SCENIC/pyscenic_mouse.sh {os.path.join(output,"sce.loom")} {os.path.join(output,"grn_result.csv")} {os.path.join(output,"ctx_result.csv")} {os.path.join(output,"aue_result.loom")}')