# Basic R Codes for Bioinformatics Analysis  
# 生信分析与绘图 R 脚本集合

This repository provides a curated collection of reusable R scripts for bioinformatics data analysis and publication-ready figure generation.

本仓库汇集了一系列可复用的 R 脚本，主要用于生信数据分析及论文级（CNS风格）高质量图形绘制。

---

## Repository Structure / 仓库结构

```
.github/workflows/   # CI workflow (auto-check & example plotting)
scripts/             # Main script collection
  ├─ ci/             # Scripts for CI testing
  ├─ enrichment/     # Enrichment analysis (e.g. GSEA)
  └─ plotting/       # Visualization scripts
data/example/        # Small example datasets
results/             # Output directory (figures/tables)
```

- **CI workflow (`figure-check.yml`)**  
  Automatically installs R dependencies, validates scripts, and renders a demo figure.  
  自动安装 R 环境、检查脚本并生成示例图。

---

## Plotting Scripts / 绘图脚本

### General visualization

- `volcano_plot.R` — Standard volcano plot  
- `volcano_yaxis_manual_or_auto.R` — Volcano plot with flexible y-axis scaling  
- `protein_bubble_single_protein.R` — Bubble plot (single protein)  
- `protein_bubble_multi_protein.R` — Bubble plot (multi-protein comparison)  
- `rank_plot_manual.R` — Ranked plot with manual labeling/highlight control  

### AML-specific analysis

- `gse116256_aml_mapping_multivolcano.R`  

  **Functionality:**
  - AML → healthy bone marrow **label transfer (Seurat)**
  - Cluster-level **cell identity mapping**
  - Core cluster (e.g. CD15-associated clusters) characterization
  - **Multi-volcano visualization (jjVolcano)** for differential expression

  **Typical use case:**
  - Identify differentiation state of AML clusters  
  - Compare AML subpopulations against normal hematopoiesis  
  - Generate publication-ready combined figures (UMAP + volcano)

---

## Enrichment Scripts / 富集分析脚本

- `run_preranked_gsea_from_deseq2.R`

  Perform preranked GSEA (Hallmark / GO:BP / Reactome) based on DESeq2 results.  
  基于 DESeq2 结果进行 preranked GSEA 分析（Hallmark / GO:BP / Reactome）。

---

## How to Use / 使用方法

1. Open the desired script in RStudio  
2. Modify input/output paths  
3. Run the script  

Outputs (if enabled) will be saved to `results/`.

1. 在 RStudio 打开脚本  
2. 修改输入/输出路径  
3. 运行脚本  

若脚本包含输出步骤，结果将保存在 `results/` 目录。

---

## Notes / 说明

- This repository is a **script collection**, not a full R package  
- 本仓库为**脚本集合**，非完整 R 包  

- Large raw data should NOT be committed  
- 不建议上传大型原始数据，仅保留示例数据  

---
