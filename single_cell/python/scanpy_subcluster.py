import gseapy as gp
import os
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import scanpy.external as sce

class recluster:
    def __init__(self,adata=None,species = None,neighbors=30,resolution=0.5,n_pcs = 10,method='leiden',batch_effect = 0,output=os.getcwd()):
        '''
        species  def is None,can use human or mouse
        neighbors 30
        resulution 0.5
        method leiden or louvain
        batch_effect 0
        output output diectory
        '''
        self.adata = adata
        self.species =species
        self.neighbors = neighbors
        self.resolution = resolution
        self.output = output
        self.method = method
        self.batch  = batch_effect
        self.n_pcs = n_pcs
    
    def sub_cluster(self):
        
        neighbors = self.neighbors
        resolution = self.resolution
        output = self.output
        method = self.method
        species = self.species
        batch_effect = self.batch
        n_pcs = self.n_pcs
        self.adata = self.adata.raw.to_adata()
        # sc.pp.normalize_total(self.adata, target_sum=1e4)
        # sc.pp.log1p(self.adata)
        self.adata.uns['log1p']['base'] = 2
        sc.pp.highly_variable_genes(self.adata,n_top_genes=2000)
        self.adata.raw = self.adata
        self.adata = self.adata[:,self.adata.var.highly_variable]
        print("neighbors:"+str(neighbors))
        print("resolution:"+str(resolution))
        print("output"+output)
        # sc.pp.log1p(self.adata)
        # sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        sc.tl.pca(self.adata, svd_solver='arpack')
        sc.pl.pca(self.adata, color='sample', show=False)
        save_fig(plt, output, file_name='pca_sample')
        sc.pl.pca_variance_ratio(self.adata, n_pcs=30, show=False)
        save_fig(plt, output, file_name='pca_variance')
        # Umap Tsne
        sc.pp.neighbors(self.adata, n_neighbors=neighbors, n_pcs=n_pcs)
        if batch_effect == 1:
            sce.pp.bbknn(self.adata, batch_key='sample', n_pcs=n_pcs)  # running bbknn 1.3.6
        if batch_effect == 2:
            sce.pp.harmony_integrate(self.adata, key ='sample',adjusted_basis="X_pca")

        sc.tl.umap(self.adata)
        if method == 'leiden':
            sc.tl.leiden(self.adata, resolution=resolution)
        else:
            sc.tl.louvain(self.adata, resolution=resolution)
        # plot cluster
        sc.pl.umap(self.adata, color=method, show=False)
        save_fig(plt, output, file_name='sub_cluster_umap')
        sc.pl.umap(self.adata, color=method,legend_loc ='on data', show=False)
        save_fig(plt, output, file_name='sub_cluster_umap_withlabel')
        sc.tl.tsne(self.adata)
        sc.pl.tsne(self.adata, color=method, size=25, show=False)
        save_fig(plt, output, file_name='sub_cluster_tsne')
        sc.pl.tsne(self.adata, color=method,size=25, legend_loc ='on data', show=False)
        save_fig(plt, output, file_name='sub_cluster_tsne_withlabel')
        # cal marker gene
        # sc.pp.log1p(self.adata)
        print(self.adata.obs.columns.tolist())
        sc.tl.rank_genes_groups(self.adata,
                                groupby=f"{method}",
                                method='wilcoxon',
                                pts=True)
        sc.pl.rank_genes_groups(self.adata, n_genes=10, sharey=False,show=False)
        save_fig(plt, output, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})
        # enrich
        find_group_marker(adata=self.adata,species=species,output=output)
        

        from matplotlib.pyplot import rc_context
        sc.set_figure_params(dpi=80, color_map = 'viridis_r')

        # group = self.adata.obs['group'].drop_duplicates().tolist()
        # if len(group) > 1:
        #     axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
        #     fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*5, 6))
        #     # keys = self.adata.obs['cell_type'].cat.categories.to_list()
        #     # values = list(self.adata.uns['cell_type_colors'])
        #     # color_dict = dict(zip(keys,values))
        #     for i in range(len(group)):
        #         if i == len(group) - 1:
        #             legend = 'right margin'
        #         else:
        #             legend = None
        #         sc.pl.umap(self.adata[self.adata.obs["group"]==group[i]], size = 20,color=method, legend_loc=legend,
        #                 #  cmap = color_dict,
        #                    title=group[i], ax=axs[i], show=False)
        #     plt.savefig(os.path.join(output, 'umap_groups.png'), facecolor='white', bbox_inches="tight", dpi=300)
        #     plt.savefig(os.path.join(output, 'umap_groups.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
        group = self.adata.obs['group'].drop_duplicates().tolist()
        if len(group) > 1:
            axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
            fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
            keys = self.adata.obs[method].cat.categories.to_list()
            values = list(self.adata.uns[method+'_colors'])
            color_dict = dict(zip(keys,values))
            for i in range(len(group)):
                tmp = self.adata[self.adata.obs["group"]==group[i]]
                sc.pl.umap(tmp, size = 20,color=method, 
                        legend_loc='on data',palette = color_dict, title=group[i], ax=axs[i], show=False)
        save_fig(plt,out_dir=output,file_name='umap_split_groups', fig_size={'w':len(group)*6, 'h':6})
        self.adata.write_h5ad(os.path.join(output,"result.h5ad"))
        return self.adata
    
    

    
# enrich 
def find_group_marker(adata, species = 'human' , output =os.getcwd() , logfc=1, fdr = 0.05, pct = 0.25):
        '''
        '''
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        all_diff = []
        for c in groups:
            out = os.path.join(output, f'cluster{c}')
            os.makedirs(out, exist_ok=True)
            tmp = sc.get.rank_genes_groups_df(adata, group=c)
            # tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
            tmp = tmp[(tmp['logfoldchanges']>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
            tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
            tmp['cluster'] = c
            all_diff.append(tmp)
            try:
                if list(tmp['names']) != []:
                    if species == 'human':
                        enr_go = gp.enrichr(gene_list=list(tmp['names']),
                            gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                            # background='hsapiens_gene_ensembl',
                            organism='Human', 
                            cutoff=0.05,
                            no_plot=False,
                            outdir=out,
                            verbose=False
                        )
                    if species == 'mouse':
                        enr_go = gp.enrichr(gene_list=list(tmp['names']),
                            gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse'],
                            # background='hsapiens_gene_ensembl',
                            organism='Mouse', 
                            cutoff=0.05,
                            no_plot=False,
                            outdir=out,
                            verbose=False
                            )
            except:
                try:
                    if list(tmp['names']) != []:
                        if species == 'human':
                            enr_go = gp.enrichr(gene_list=list(tmp['names']),
                                gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                                # background='hsapiens_gene_ensembl',
                                organism='Human', 
                                cutoff=0.05,
                                no_plot=True,
                                outdir=out,
                                verbose=False
                            )
                        if species == 'mouse':
                            enr_go = gp.enrichr(gene_list=list(tmp['names']),
                                gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse'],
                                # background='hsapiens_gene_ensembl',
                                organism='Mouse', 
                                cutoff=0.05,
                                no_plot=True,
                                outdir=out,
                                verbose=False
                                )
                except:
                    continue

        all_diff = pd.concat(all_diff, axis=0)
        all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        all_diff.to_csv(os.path.join(output, 'all_cluster_marker.xls'), sep='\t', index=None)
        all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(output, 'cluster_diff_statistics.xls'), sep='\t')

def draw_cell_type(adata,groupby='cell_type',markerlist = [],outdir = os.getcwd()):
    '''
    adata : input anndata
    groupby : celltype name
    markerlist : list of cell type marker gene
    '''
    sc.pl.umap(adata,color= groupby,size=20,show=False)
    save_fig(plt,out_dir=outdir,file_name='celltype_umap',fig_size={'w':7, 'h':5})
    sc.pl.umap(adata,color= groupby,size=20,legend_loc='on data',show=False)
    save_fig(plt,out_dir=outdir,file_name='celltype_umap_withlabel',fig_size={'w':5, 'h':5})
    sc.pl.dotplot(adata,markerlist,groupby=groupby,show=False)
    save_fig(plt,out_dir=outdir,file_name='celltype_dotplot',fig_size={'w':len(markerlist)/4*5, 'h':5})
    group = adata.obs['group'].drop_duplicates().tolist()
    if len(group) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
        fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
        keys = adata.obs[groupby].cat.categories.to_list()
        values = list(adata.uns[groupby+'_colors'])
        color_dict = dict(zip(keys,values))
        for i in range(len(group)):
            tmp = adata[adata.obs["group"]==group[i]]
            sc.pl.umap(tmp, size = 20,color=groupby, 
                       legend_loc='on data',palette = color_dict, title=group[i], ax=axs[i], show=False)
    save_fig(plt,out_dir=outdir,file_name=f'{groupby}_umap_split_groups', fig_size={'w':len(group)*6, 'h':6})
    adata.obs.to_csv(os.path.join(outdir,'cell_type_meta.xls'),sep="\t")
    adata.write_h5ad(os.path.join(outdir,'Result_celltype.h5ad'))






def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()