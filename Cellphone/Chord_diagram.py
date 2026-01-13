import os
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
from cellphonedb.utils import search_utils
os.chdir('/lvdata/wzb/scRNA/FW2023-656/2024.01.19/Cellphone/WT')

adata = ad.read_h5ad('/lvdata/wzb/scRNA/FW2023-656/2024.01.19/wt_Result_human.h5ad')

means = pd.read_csv('statistical_analysis_means_.txt',sep="\t")
pvals = pd.read_csv('statistical_analysis_pvalues_.txt',sep="\t")
decon = pd.read_csv('statistical_analysis_deconvoluted_.txt',sep="\t")
significant_means = pd.read_csv('statistical_analysis_significant_means_.txt',sep="\t")

cell_type_list = adata.obs['cell_type'].drop_duplicates().to_list()
for i in cell_type_list:

    search_utils_result = search_utils.search_analysis_results(
        query_cell_types_1 =[i],
        query_cell_types_2 = cell_type_list,
        significant_means = significant_means,
        deconvoluted = decon,
        separator = '|',                                                # separator (default: |) employed to split cells (cellA|cellB).
        long_format = False   
    )
    search_utils_result['Total'] = search_utils_result.iloc[:, 5:].sum(axis=1)
    df_sorted = search_utils_result.sort_values(by='Total', ascending=False)
    top_3_rows = df_sorted.head(3)
    # gene_columns = ['gene_a', 'gene_b']
    gene_columns = ['gene_a']
    top_3_genes = top_3_rows[gene_columns].values.flatten().tolist()
    top_3_genes_unique = list(set(top_3_genes))
    top_3_genes_unique =[gene for gene in top_3_genes_unique if pd.notna(gene)]
    print(top_3_genes_unique)

    # p = kpy.plot_cpdb_chord(adata,means,pvals,decon,celltype_key='cell_type',cell_type1=i,cell_type2='.',genes=['CNR1'])
    p = kpy.plot_cpdb_chord(adata,means,pvals,decon,celltype_key='cell_type',cell_type1=i,cell_type2='.',
                            remove_self =True,genes=top_3_genes_unique,figsize=(6, 6),gap=0,
                            labelposition=50)
    p.save(f"{i}_chord_diagram",format="png",dpi=300)
    p.save(f"{i}_chord_diagram",dpi=300)

