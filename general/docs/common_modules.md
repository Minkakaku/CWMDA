# 通用模块说明

本目录用于存放所有组学可复用的公共模块，避免重复开发。

## 模块列表

- `omics/general/scripts/config_utils.py`
  - `load_config(common_path, specific_path)`：读取通用配置 + 组学配置并合并。
  - `resolve_path(value, base_dir=None)`：解析相对路径。
- `omics/general/scripts/log_utils.py`
  - `setup_logger(name, log_dir, level)`：统一日志输出到控制台 + 文件。
- `omics/general/scripts/file_utils.py`
  - `ensure_dir(path)`：创建输出目录。
  - `validate_file(path, label)`：检查输入文件是否存在。
- `omics/general/scripts/qc_common.py`
  - `write_qc_report(stats, output_dir, tsv_name, html_name)`：统一输出 TSV + HTML 报告。

## 使用示例

```python
from omics.general.scripts.config_utils import load_config
from omics.general.scripts.qc_common import write_qc_report

config = load_config(
    "omics/general/config/common_config.yaml",
    "omics/genomics/config/genome_config.yaml",
)
write_qc_report({"sample": "demo"}, "data/genome/output", "qc_stats.tsv", "qc_report.html")
```
