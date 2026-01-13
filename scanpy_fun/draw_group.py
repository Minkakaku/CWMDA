import os
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()
def draw_group(adata,out_dir = os.getcwd(),method = 'leiden',size=20):
    '''
    '''
    # group = adata.obs['group'].drop_duplicates().tolist()
    group = adata.obs['group'].cat.categories.to_list()
    if len(group) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
        fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
        keys = adata.obs[method].cat.categories.to_list()
        values = list(adata.uns[method+'_colors'])
        color_dict = dict(zip(keys,values))
        for i in range(len(group)):
            tmp = adata[adata.obs["group"]==group[i]]
            sc.pl.umap(tmp, size = 20,color=method, 
                       legend_loc='on data',palette = color_dict, title=group[i], ax=axs[i], show=False)
    save_fig(plt,out_dir=out_dir,file_name='umap_split_groups', fig_size={'w':len(group)*6, 'h':6})
    if len(group) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
        fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
        keys = adata.obs[method].cat.categories.to_list()
        values = list(adata.uns[method+'_colors'])
        color_dict = dict(zip(keys,values))
        for i in range(len(group)-1):
            tmp = adata[adata.obs["group"]==group[i]]
            sc.pl.umap(tmp, size = size,color=method, 
                        legend_loc='none',palette = color_dict, title=group[i], ax=axs[i], show=False)
        i  = len(group)-1
        tmp = adata[adata.obs["group"]==group[i]]
        sc.pl.umap(tmp, size = size,color=method, 
                    legend_loc='right margin',palette = color_dict, title=group[i], ax=axs[i], show=False)
    save_fig(plt,out_dir=out_dir,file_name='umap_split_groups_withoutlegend', fig_size={'w':len(group)*6, 'h':6})
    # if len(group) > 1:
    #     axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
    #     fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
    #     keys = adata.obs[method].cat.categories.to_list()
    #     values = list(adata.uns[method+'_colors'])
    #     color_dict = dict(zip(keys,values))
    #     for i in range(len(group)):
    #         tmp = adata[adata.obs["group"]==group[i]]
    #         sc.pl.tsne(tmp, size = 20,color=method, 
    #                    legend_loc='on data',palette = color_dict, title=group[i], ax=axs[i], show=False)
    # save_fig(plt,out_dir=output,file_name='tsne_split_groups', fig_size={'w':len(group)*6, 'h':6})

def getcolor(adata,output = os.getcwd(),method = 'leiden'):

    keys = adata.obs[method].cat.categories.to_list()
    values = list(adata.uns[method+'_colors'])
    color_dict = dict(zip(keys,values))
    df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
    df.columns.names = ["cell_type"]
    df.to_csv(os.path.join(output,method+'color_dict.xls'),sep="\t",index=True)