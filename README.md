# Basic R Codes for Bioinformatics Analysis / 生信分析与绘图 R 脚本集合

This repository collects reusable R scripts for bioinformatics data analysis and publication-ready figures.  
本仓库收集可复用的 R 脚本，侧重生信数据分析与用于论文发表（CNS风格）的高质量绘图复现。

## Repository structure / 仓库结构

- `.github/workflows/figure-check.yml`  
  GitHub Actions workflow: install R & plotting packages, check scripts, render a CI example figure.  
  GitHub Actions 工作流：自动安装 R/绘图包、检查脚本并渲染一个 CI 示例图。

- `scripts/ci/`  
  Scripts specifically for CI checks (example figure rendering).  
  专用于 CI 检查的脚本（示例图渲染）。

- `scripts/enrichment/`  
  Enrichment analysis scripts (e.g., preranked GSEA from DESeq2 results).  
  富集分析脚本（如基于 DESeq2 的 preranked GSEA）。

- `scripts/plotting/`  
  Plotting scripts for commonly used figures.  
  常用绘图脚本。

- `data/example/`  
  Small example datasets used by scripts.  
  脚本用的小示例数据。

- `results/`  
  Output directory for generated figures/tables (not all outputs are committed).  
  输出目录（生成的图/表；不一定全部提交）。

## Plotting scripts / 绘图脚本

- `scripts/plotting/volcano_plot.R` — Volcano plot  
- `scripts/plotting/volcano_yaxis_manual_or_auto.R` — Volcano plot with manual/auto y-axis options  
- `scripts/plotting/protein_bubble_single_protein.R` — Protein bubble plot (single protein)  
- `scripts/plotting/protein_bubble_multi_protein.R` — Protein bubble plot (multi-protein)  
- `scripts/plotting/rank_plot_manual.R` — Ranked plot with manual highlight/label logic  

- `scripts/plotting/gse116256_aml_mapping_multivolcano.R`
— GSE116256 AML cluster-to-healthy mapping + core-cluster multi-volcano plots（Seurat Label Transfer + jjVolcano）
## Enrichment scripts / 富集分析脚本

- `scripts/enrichment/run_preranked_gsea_from_deseq2.R`  
  Run Hallmark / GO:BP / Reactome GSEA based on DESeq2 output.  
  基于 DESeq2 输出，跑 Hallmark / GO:BP / Reactome GSEA。

## How to use / 使用方法

1. Open a script in RStudio and adjust input/output paths as needed.  
2. Run the script.  
3. Generated outputs will be saved in `results/` (if the script writes outputs).  

1. 在 RStudio 打开脚本，根据需要调整输入/输出路径。  
2. 执行脚本。  
3. 若脚本保存输出，输出会写入 `results/` 目录。  

## Notes / 说明

- This repo is a script collection rather than a full R package.  
- 本仓库是脚本集合，而不是完整 R package。  

- Do not commit large raw data files; keep only small example files.  
- 不要提交超大的原始数据；只保留小示例数据即可。  
