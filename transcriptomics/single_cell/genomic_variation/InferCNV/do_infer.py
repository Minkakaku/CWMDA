import pandas as pd
import infercnvpy as cnv
import os
import matplotlib.pyplot as plt
import scanpy as sc
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

gene_pos = pd.read_csv('/lvdata/wzb/pipline/InferCNV/inferCNV_R/mm10_2020A.txt',sep="\t",index_col=0,header=None)
gene_pos.columns = ["chromosome", "start", "end"]
gene_pos.index.name = None
for i in ["chromosome", "start", "end"]:
    adata.var[i] = gene_pos[i]
    
cnv.tl.infercnv(
    adata,
    reference_key='cell_type',
    reference_cat= ['B cells',  'Myeloid cells', 'Neutrophils', 'T cells'],
    window_size=250
)

outdir = '/lvdata/wzb/scRNA/FWSC20240564/2024.10.08/inferCNV'
os.makedirs(outdir,exist_ok=True)
cnv.pl.chromosome_heatmap(adata, groupby="cell_type",show=False)
save_fig(plt,out_dir=outdir,file_name='chromosome_heatmap',fig_size={'w':10,'h':7})
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
adata.obs['cellname'] = adata.obs.index
cnv.tl.cnv_score(adata,groupby='cellname')
# cnv.tl.cnv_score(adata)

sc.pl.umap(
    adata,
    color=["cnv_score",'cell_type','cnv_leiden','group'],
    show=False,
    size= 20,
    wspace = 0.5
)
save_fig(plt,out_dir=outdir,file_name='CNV_score',fig_size={'w':23,'h':7})
sc.pl.dotplot(adata,'cnv_score',groupby='cnv_leiden')
sc.pl.umap(adata,color='cnv_leiden',legend_loc='on data')