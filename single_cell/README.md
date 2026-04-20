# Single-Cell Codebase

这里是仓库内单细胞代码的统一新目录。重构后的目标是：

- 按功能拆分，而不是按历史项目拆分
- 高频脚本统一成可传参入口
- 小工具尽量函数化，避免写死项目路径
- 旧目录迁移后，只保留 `single_cell` 作为单细胞代码入口

## 目录结构

- `python`
  - `run_scanpy_pipeline.py` / `scanpy_pipeline.py`：Scanpy 主流程
  - `sample_discovery.py`：10x 样本自动发现
  - `infercnv_pipeline.py`：InferCNV Python 入口
  - `cellphonedb_pipeline.py`：CellPhoneDB 入口
  - `scenic_prepare.py`：SCENIC 前处理
  - `scanpy_draw_group.py` / `scanpy_find_marker.py` / `scanpy_subcluster.py`：Scanpy 小工具
- `r`
  - `seurat_main_*`：Seurat 主流程
  - `seurat_subcluster_*`：Seurat 子群流程
  - `cellchat_pipeline.R`：CellChat
  - `monocle2_*` / `monocle3_pipeline.R` / `slingshot_pipeline.R` / `cytotrace_pipeline.R`：轨迹与发育分析
  - `aucell_pipeline.R` / `nichenetr_pipeline.R` / `scmetabolism_pipeline.R` / `enrichment_pipeline.R` / `louper_export.R`：功能分析与导出
  - `proportion_pipeline.R` / `roe_pipeline.R`：组成比例与 Ro/e
  - `common_utils.R`：R 侧公共参数解析、输出保存和目录工具
- `sh`
  - `cellranger_gex.sh`
  - `pyscenic_human.sh`
  - `pyscenic_mouse.sh`
  - `velocyto_run10x.sh`
  - `velocyto_cleanup.sh`
  - `setup_seurat_env.sh`

## 典型入口

### 1. CellChat

```bash
Rscript single_cell/r/cellchat_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/cellchat \
  --group-col group \
  --celltype-col cell_type \
  --species human
```

### 2. NicheNetR

```bash
Rscript single_cell/r/nichenetr_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/nichenetr \
  --celltype-col cell_type \
  --condition-col group \
  --receiver Fibroblast \
  --sender-celltypes Myeloid,Tcell \
  --condition-oi Case \
  --condition-reference Control \
  --resources-dir /path/to/nichenetr_resources
```

### 3. 比例和 Ro/e

```bash
Rscript single_cell/r/proportion_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/proportion \
  --celltype-col cell_type \
  --group-col group

Rscript single_cell/r/roe_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/roe \
  --celltype-col cell_type \
  --group-col group
```

### 4. Monocle3

```bash
Rscript single_cell/r/monocle3_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/monocle3 \
  --cluster-col cell_type \
  --root-col cell_type \
  --root-value Naive
```

### 5. AUCell

```bash
Rscript single_cell/r/aucell_pipeline.R \
  --input-rds your_seurat.rds \
  --out-dir result/aucell \
  --gene-set-file gene_sets.tsv \
  --ident-col cell_type \
  --split-col group
```

## 参数风格

重构后的脚本优先使用统一的 `--key value` 风格，而不是旧的固定位置参数。常见参数包括：

- `--input-rds`
- `--out-dir`
- `--celltype-col`
- `--group-col`
- `--species`
- `--assay`

## 说明

- 这次重构以“统一入口、去硬编码、保留算法逻辑”为主。
- Monocle2 相关脚本保留了原有主体逻辑，但去掉了固定环境路径。
- 当前环境里没有可直接调用的 Python 解释器，所以 Python 侧 CLI 只能做静态整理，没法在本机逐个实际运行验证。
