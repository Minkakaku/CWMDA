import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import cosg as cosg
import numpy as np
import gseapy as gp
import scanpy as sc
import os
import re
def find_marker_gene(adata, group,species=None,out_dir=os.getcwd(),filter = True, logfoldchanges=0.25,pvals_adj=0.05,pct_nz_group=0.1, method='wilcoxon', reference = 'rest',gene_annot=None):
    """
    查找标记基因。

    Args:
        adata: 输入的 AnnData 对象。
        group: 细胞类型群组（例如 leiden）。
        species (str, optional): 物种，如 'human' 或 'mouse'。默认为 None。
        out_dir (str, optional): 输出目录。默认为当前工作目录。
        filter (bool, optional): 是否过滤。默认为 True。
        logfoldchanges (float, optional): 对数折叠变化阈值。默认为 0.25。
        pvals_adj (float, optional): 调整后的 p 值阈值。默认为 0.05。
        pct_nz_group (float, optional): 组中非零百分比阈值。默认为 0.1。
        method (str, optional): 使用的方法，如 'wilcoxon'。默认为 'wilcoxon'。
        reference (str, optional): 参考组。默认为 'rest'。
        gene_annot (optional): 基因注释。默认为 None。

    Returns:
        写出group类的所有marker gene。
    """
    sample_counts = adata.obs[group].value_counts()
    problematic_groups = sample_counts[sample_counts == 1].index.to_list()
    adata = adata[~adata.obs[group].isin(problematic_groups)]
    adata.uns['log1p']['base'] = None
    sc.tl.rank_genes_groups(adata, group, method=method,reference= reference,pts=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    all_diff = []
    for c in groups:
        
        # c = 'CD4+TN'
        d = re.sub(r'[+/\\-]',"_",c)
        out = os.path.join(out_dir, f'cluster{d}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        if filter :
            tmp = tmp[(tmp['logfoldchanges']>logfoldchanges) & (tmp['pvals_adj']<pvals_adj) & (tmp['pct_nz_group'] >=pct_nz_group)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)
        try:
            if list(tmp['names']) != []:
                if species == 'human':
                    enr_go = gp.enrichr(gene_list=list(tmp['names']),
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human','Reactome_2016'],
                        # background='hsapiens_gene_ensembl',
                        organism='Human', 
                        cutoff=0.05,
                        no_plot=False,
                        outdir=out,
                        verbose=False
                    )
                if species == 'mouse':
                    enr_go = gp.enrichr(gene_list=list(tmp['names']),
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse','Reactome_2016'],
                        # background='mmusculus_gene_ensembl',
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
                            # ackground='mmusculus_gene_ensembl',
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
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')
    # sc.pl.rank_genes_groups_dotplot(adata,groups=group,n_genes=5,groupby=group,show=False,save=f'{out_dir}/bubbleplot.pdf')