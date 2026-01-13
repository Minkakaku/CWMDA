import os
import anndata
import argparse
import warnings
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(prog='pipeline_cellphonedb.py', description="cellphonedb 分析流程")
parser.add_argument('-i', action='store', dest='h5ad', required=True, help='adata.h5ad')
parser.add_argument('-o', action='store', dest='out_dir', help='Result output directory. default:cellphonedb_result')
parse = parser.parse_args()

if not parse.out_dir:
    parse.out_dir = 'cellphonedb_result'
if os.path.exists(parse.out_dir):
    os.system('rm -r {0}'.format(parse.out_dir))
os.makedirs(parse.out_dir, exist_ok=True)

cpdb_file_path = '/home/hf/project/shen/mfcellphonedb/v5.0.0/cellphonedb.zip'
counts_file_path = parse.h5ad
out_path = parse.out_dir

adata = anndata.read_h5ad(counts_file_path)
# gene = adata.var
# gene.to_csv(os.path.join(out_path,'gene.tsv'), index=True, sep = '\t')
meta = adata.obs[['cell_type']]
meta.to_csv(os.path.join(out_path,'meta.tsv'), index=True, sep = '\t')
meta_file_path = os.path.join(out_path,'meta.tsv')

#adata.X.to_csv(os.path.join(out_path,'adata.tsv'),index=True, sep = '\t')
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    #microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis. 1000
    threshold = 0.01,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 4,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 15000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = ""                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )


kpy.plot_cpdb_heatmap(
        adata = adata,
        pvals = pvalues,
        celltype_key = "cell_type",
        figsize = (15,15),
        title = "Sum of significant interactions",
        symmetrical = False,
    )
plt.savefig(os.path.join(out_path,'heatmap.png'),facecolor='white', bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_path,'heatmap.pdf'),facecolor='white')

cell_type = adata.obs['cell_type'].drop_duplicates().to_list()
cell_type = adata.obs['cell_type'].drop_duplicates().to_list()
from cellphonedb.utils import search_utils

# cell_type
for i in cell_type:
    try:
        search_utils_result = search_utils.search_analysis_results(
            query_cell_types_1 =[i],
            query_cell_types_2 = cell_type,
            significant_means = significant_means,
            deconvoluted = deconvoluted,
            separator = '|',                                                # separator (default: |) employed to split cells (cellA|cellB).
            long_format = False   
        )
        search_utils_result['Total'] = search_utils_result.iloc[:, 5:].sum(axis=1)
        df_sorted = search_utils_result.sort_values(by='Total', ascending=False)
        top_10_rows = df_sorted.head(10)
        gene_columns = ['gene_a', 'gene_b']
        top_10_genes = top_10_rows[gene_columns].values.flatten().tolist()
        top_10_genes_unique = list(set(top_10_genes))
        top_10_genes_unique =[gene for gene in top_10_genes_unique if pd.notna(gene)]
        print(top_10_genes_unique)
    except:
        top_10_genes = None
    if top_10_genes != None:
        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = i,
            cell_type2 = ".", 
            means = means,
            pvals = pvalues,
            genes = top_10_genes_unique,
            celltype_key = "cell_type",             #细胞类型的列名称
            figsize = (20,10),
            title = "Interactions",
            keep_significant_only = True,
            max_size = 6,
            highlight_size = 0.75,
            standard_scale = True
        )

        p.save(os.path.join(out_path,f'{i}_Top10_bubble.png'),facecolor='white', bbox_inches='tight', dpi=300)
        p.save(os.path.join(out_path,f'{i}_Top10_bubble.pdf'),facecolor='white')
    else:
        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = i,
            cell_type2 = ".", 
            means = means,
            pvals = pvalues,
            # genes = top_10_genes_unique,
            celltype_key = "cell_type",             #细胞类型的列名称
            figsize = (20,10),
            title = "Interactions",
            keep_significant_only = True,
            max_size = 6,
            highlight_size = 0.75,
            standard_scale = True
        )

        p.save(os.path.join(out_path,f'{i}_bubble.png'),facecolor='white', bbox_inches='tight', dpi=300)
        p.save(os.path.join(out_path,f'{i}_bubble.pdf'),facecolor='white')

os.system('rm  {0}'.format(meta_file_path))



