from operator import index
from re import T
import streamlit as st
import scanpy as sc
import base64
import os
import random
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


data_list = os.listdir('/home/lst/scRNA_streamlit/cell_cell_Communication/result')
project_list=st.multiselect("projects",data_list,default=['test3'])
merge_cc,all_cell_types_list= get_cell_pair(project_list)
all_cell_types_list=st.multiselect("interacting_cell type",all_cell_types_list)
project_name= st.text_input('请输入项目名：')

submitted = st.button("🚀提交")
all_cell_types_list = []
if submitted:
    if project_name:
        out_dir = os.path.join('/home/lst/cell_cell_Communication/result/CellPhoneDB_group', project_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            st.success(f'已成功创建文件夹：{out_dir}')
            out_dir = os.path.join('/home/lst/cell_cell_Communication/result/CellPhoneDB_group', project_name)
        else:
            random_chars = '_'.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 4)) 
            new_title = project_name+random_chars
            out_dir = os.path.join('/home/lst/cell_cell_Communication/result/CellPhoneDB_group', new_title)
            st.info('该项目名已存在，新项目保存为'+new_title, icon="ℹ️")
    filtered_df = merge_cc[merge_cc['variable'].apply(lambda x: any(cell_type in x for cell_type in all_cell_types_list))]
    filtered_df.to_csv(os.path.join(out_dir,'Mean_pvalue.txt'),sep='\t',index=False)
    os.system('Rscript ./bubble.R -i {0} -o {1}'.format(os.path.join(out_dir,'Mean_pvalue.txt'),out_dir))

