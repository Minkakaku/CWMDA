#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from operator import index
from re import T
import streamlit as st
import scanpy as sc
import base64
import os
import pandas as pd
import scanpy as sc
import ktplotspy as kpy
from pagess.plot_cpdb import plot_cpdb
import matplotlib.pyplot as plt
from PIL import Image
import cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from pagess.functions import cell_select,get_mean_p,get_cell_pair
for k, v in st.session_state.items():
    st.session_state[k] = v
st.session_state.items()

if st.session_state['_state_']['adata'] is not None:
    st.markdown('''##### 分析参数''')
    adata=st.session_state['_state_']['adata']

    out_dir = os.path.join(st.session_state['out_dir'], 'CellPhoneDB')
    os.makedirs(out_dir, exist_ok=True)
    obs_list = list(adata.obs['cell_type'].unique())
    obs_list = [item for item in obs_list if isinstance(item, str)]
    st.info('选择细胞类型', icon="ℹ️")
    cell_type_list=st.multiselect("Cell type", obs_list,default=obs_list)
    #adata_cpdb = adata[adata.obs['cell_type'].isin(cell_type_list), :]
    adata_cpdb = cell_select(adata,cell_type_list,'cpdb',out_dir)
    st.session_state['_state_']['adata_cpdb']=adata_cpdb
    
    meta = adata_cpdb.obs[['cell_type']]
    meta.to_csv(os.path.join(out_dir,'meta.tsv'), index=True, sep = '\t')
    counts_file_path = os.path.join(out_dir,'adata_cpdb.h5ad')
    meta_file_path = os.path.join(out_dir,'meta.tsv')

    col1, col2, col3 = st.columns(3)
    st.session_state['_state_']['min_pct_cell'] = col1.number_input(label='过滤低于n%细胞表达的基因：', min_value=0, max_value=100, value=st.session_state['_state_']['min_pct_cell'], step=1)
    st.session_state['_state_']['cpdb_p'] = col2.number_input(label='Pvalue：', min_value=0.0, max_value=1.0, value=st.session_state['_state_']['cpdb_p'], step=0.01)
    st.session_state['_state_']['subsampling'] = st.checkbox('子采样', value=st.session_state['_state_']['subsampling'], help='子采样，当细胞数较多时建议进行子采样')
    st.session_state['_state_']['subsampling_log'] = st.checkbox('子采样对数转换', value=st.session_state['_state_']['subsampling_log'], help='当输入数据未进行对数转换时，进行子采样勾选对数转换')
    st.session_state['_state_']['subsampling_num_pc'] = col3.number_input(label='子采样PC数', min_value=0, max_value=200, value=st.session_state['_state_']['subsampling_num_pc'], step=1)
    #st.session_state['_state_']['symmetrical'] = st.checkbox('热图是否对称显示', value=st.session_state['symmetrical'], help='不区分配受体')
    submitted = st.button("🚀提交")
    if submitted:

        deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
            cpdb_file_path = './db/cellphonedb.zip',                        # mandatory: CellPhoneDB database zip file.
            meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
            counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
            counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
            #microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
            iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
            threshold = st.session_state['_state_']['min_pct_cell']/100,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
            threads = 4,                                     # number of threads to use in the analysis.
            debug_seed = 42,                                 # debug randome seed. To disable >=0.
            result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
            pvalue = st.session_state['_state_']['cpdb_p'],                                   # P-value threshold to employ for significance.
            subsampling = st.session_state['_state_']['subsampling'],                             # To enable subsampling the data (geometri sketching).
            subsampling_log = st.session_state['_state_']['subsampling_log'],                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
            subsampling_num_pc = st.session_state['_state_']['subsampling_num_pc'],                        # Number of componets to subsample via geometric skectching (dafault: 100).
            #subsampling_num_cells = 3000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
            separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
            debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
            output_path = out_dir,                          # Path to save results.
            output_suffix = 'result'                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
            )


        st.session_state['_state_']['symmetrical'] = st.checkbox('热图是否对称显示', value=st.session_state['_state_']['symmetrical'], help='不区分配受体')
        kpy.plot_cpdb_heatmap(
            adata = adata_cpdb,
            pvals = pvalues,
            log1p_transform= False,
            celltype_key = 'cell_type',
            figsize = (6,6),
            title = "Sum of significant interactions",
            symmetrical = st.session_state['_state_']['symmetrical']
        )
        plt.savefig(os.path.join(out_dir,'heatmap.png'),facecolor='white', bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(out_dir,'heatmap.pdf'),facecolor='white')


        cell_type1=st.multiselect("Cell type1", obs_list,default=obs_list)
        cell_type2=st.multiselect("Cell type2", obs_list,default=obs_list)
        p = plot_cpdb(
        adata = adata_cpdb,
        cell_type1 = '.',
        cell_type2 = '.', 
        means = means,
        pvals = pvalues,
        celltype_key = 'cell_type',
        figsize = (20,15),
        title = "Interactions",
        max_size = 4,
        #highlight_size = 0.3,
        standard_scale = True
        )

        p.save(os.path.join(out_dir,'bubble.png'),facecolor='white', bbox_inches='tight', dpi=300)
        p.save(os.path.join(out_dir,'bubble.pdf'),facecolor='white')
        heatmap = Image.open(os.path.join(out_dir,'heatmap.png'))
        bubble = Image.open(os.path.join(out_dir,'bubble.png'))
        tabs = st.tabs(["🚀Heatmap plot", "📈bubble"])
        tabs[0].image(heatmap, caption='', use_column_width=True)
        tabs[1].image(bubble, caption='', use_column_width=True)
        df_p_m = get_mean_p(out_dir)
        df_p_m.to_csv(os.path.join(out_dir,'Mean_pvalue.txt'),sep='\t',index=False)
    st.markdown('''
    ##### 请先上传数据
    ''')