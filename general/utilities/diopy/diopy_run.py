import anndata as ad
import diopy
import argparse
import warnings
import os

warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(prog='diopy.py', description="h5ad2h5")
parser.add_argument('-i', action='store', dest='h5ad', required=True, help='adata.h5ad')
parser.add_argument('-o', action='store', dest='out_dir',required=True, help='Result output directory.')
parser.add_argument('-n', action='store', dest='name',required=True, help='Result output h5name.')
parser = parser.parse_args()


adata = ad.read_h5ad(parser.h5ad)
os.chdir(parser.out_dir)
adata = adata.raw.to_adata()
adata.raw = adata
adata.obs_names_make_unique()
# 这里有时候会报错
# diopy.output.write_h5(adata,file = parser.name,save_X=True,save_graph = True)
diopy.output.write_h5(adata,file = parser.name,save_X=False,save_graph = True)
# print(f"{os.path.join(parser.out_dir,parser.name)}")
os.system(f'sh /data/NAS03/backup/wzb/wzb/pipline/dio/docker_run.sh "{os.path.join(parser.out_dir,parser.name)}" "{parser.out_dir}" "h5_trans.rds"')