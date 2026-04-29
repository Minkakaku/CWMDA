# Python Single Cell Workflow

这里放 Python/Scanpy 版本的单细胞主流程代码。文件已经按 R 主流程的结构整理为 `main/`、`subcluster/` 和 `tools/`，同时根目录保留了兼容入口，原来的主要运行方式不需要改。

## 目录

- `main/`
  - `run_scanpy_pipeline.py`：Scanpy 主流程入口。
  - `scanpy_pipeline.py`：主流程函数。
  - `scanpy_find_marker.py`：marker 工具。
  - `scanpy_draw_group.py`：分组绘图工具。
- `subcluster/`
  - `scanpy_subcluster.py`：Scanpy 亚群分析。
- `tools/`
  - `sample_discovery.py`：自动发现 10x 样本并按样本名推断分组。
  - `io_utils.py`：通用读写工具。
  - `mouse_to_human.py`：鼠/人基因名转换工具。
- 根目录同名 `.py` 文件是兼容入口，会转发到整理后的真实脚本。

## 主流程运行

原来的入口仍然可用：

```bash
python "single cell/Python single cell/run_scanpy_pipeline.py" all \
  --sample-sheet sample_sheet.tsv \
  --out-root scanpy_result
```

也可以直接调用整理后的脚本：

```bash
python "single cell/Python single cell/main/run_scanpy_pipeline.py" all \
  --sample-sheet sample_sheet.tsv \
  --out-root scanpy_result
```

## 样本自动识别

Python 版保留原命令；同等功能也已经迁移到 R：`single cell/R single cell/tools/sample_discovery.R`。

```bash
python "single cell/Python single cell/sample_discovery.py" discover \
  --data-root data2analysis \
  --output-manifest sample_manifest.tsv

python "single cell/Python single cell/sample_discovery.py" suggest \
  --manifest sample_manifest.tsv \
  --output-report group_suggestions.tsv

python "single cell/Python single cell/sample_discovery.py" assign \
  --manifest sample_manifest.tsv \
  --group-count 4 \
  --output-sample-sheet sample_sheet.tsv \
  --output-report group_suggestions.tsv
```

## 非主流程工具

CellPhoneDB、inferCNV、SCENIC 等工具已放到 `single cell/alternative tools/Python/`。根目录保留兼容入口，例如：

```bash
python "single cell/Python single cell/cellphonedb_pipeline.py" --help
python "single cell/Python single cell/infercnv_pipeline.py" --help
python "single cell/Python single cell/scenic_prepare.py" --help
```
