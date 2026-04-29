# Alternative Tools

这里存放主流程之外的单细胞和空间分析工具，例如拟时序、细胞通讯、功能打分、CNV、空间转录组、导出和上游辅助脚本。

## R

- `R/trajectory/`：Monocle2、Monocle3、Slingshot、CytoTRACE 等拟时序和发育相关分析。
- `R/communication/`：CellChat、NicheNet。
- `R/function_analysis/`：AUCell、scMetabolism、富集分析。
- `R/proportion/`：比例和 Ro/e 的旧版独立工具。
- `R/cnv/`：inferCNV R 流程和说明。
- `R/export/`：h5ad、Loupe 等导出工具。
- `R/legacy_seurat/`：旧版 Seurat 主流程和亚群流程脚本，作为参考保留。
- `R/utils/`：旧版工具函数。

## Python

- `Python/communication/`：CellPhoneDB 和弦图。
- `Python/cnv/`：inferCNV Python 相关脚本。
- `Python/scenic/`：SCENIC 前处理。

## Shell And Spatial

- `shell/workflow_scripts/`：cellranger、pyscenic、velocyto 相关 shell 脚本。
- `spatial/scripts/`：Xenium/空间转录组分析脚本。
