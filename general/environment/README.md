# 依赖环境说明

## 通用依赖

```bash
conda env create -f omics/general/environment/common_env.yml
```

## 基因组专属依赖

```bash
conda env create -f omics/genomics/environment/genome_env.yml
```

## 转录组专属依赖

```bash
conda env create -f omics/transcriptomics/environment/transcriptome_env.yml
```

> 若只运行某一组学流程，只需安装通用依赖 + 对应组学专属依赖即可。
