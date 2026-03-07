# =========================================================
# Run common preranked GSEA from DESeq2 results
# 基于DESeq2结果一次跑常用的 preranked GSEA
# =========================================================

## -------------------------
## 0. Install packages if needed
## 0. 如需安装包
## -------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_cran <- c("readr", "dplyr", "msigdbr")
pkgs_bioc <- c("clusterProfiler", "enrichplot")

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

library(readr)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

## =========================================================
## 1. Main function
## 1. 主函数
## =========================================================
run_common_gsea <- function(
    input_file,
    gene_col = "gene",
    stat_col = "stat",
    species = "Homo sapiens",
    
    run_hallmark = TRUE,
    run_gobp = TRUE,
    run_reactome = TRUE,
    
    pvalue_cutoff = 0.05,
    p_adjust_method = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    
    save_csv = TRUE,
    output_prefix = "GSEA"
) {
  
  ## -------------------------
  ## 1.1 Read input table
  ## 1.1 读取输入表
  ## -------------------------
  df <- read_csv(input_file, show_col_types = FALSE)
  
  ## Check columns
  ## 检查列名
  required_cols <- c(gene_col, stat_col)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns: ",
        paste(missing_cols, collapse = ", "),
        "\nPlease check gene_col and stat_col."
      )
    )
  }
  
  ## -------------------------
  ## 1.2 Build ranked gene list
  ## 1.2 构建排序向量
  ## -------------------------
  gene_df <- df %>%
    select(
      gene = all_of(gene_col),
      stat = all_of(stat_col)
    ) %>%
    filter(!is.na(gene), !is.na(stat)) %>%
    distinct(gene, .keep_all = TRUE)
  
  gene_list <- gene_df$stat
  names(gene_list) <- gene_df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat("Gene list prepared.\n")
  cat("Number of ranked genes:", length(gene_list), "\n\n")
  
  ## -------------------------
  ## 1.3 Prepare gene set collections
  ## 1.3 准备基因集
  ## -------------------------
  collections <- list()
  
  if (run_hallmark) {
    hallmark_df <- msigdbr(
      species = species,
      category = "H"
    ) %>%
      select(gs_name, gene_symbol)
    
    collections$Hallmark <- hallmark_df
  }
  
  if (run_gobp) {
    gobp_df <- msigdbr(
      species = species,
      category = "C5",
      subcategory = "GO:BP"
    ) %>%
      select(gs_name, gene_symbol)
    
    collections$GO_BP <- gobp_df
  }
  
  if (run_reactome) {
    reactome_df <- msigdbr(
      species = species,
      category = "C2",
      subcategory = "CP:REACTOME"
    ) %>%
      select(gs_name, gene_symbol)
    
    collections$Reactome <- reactome_df
  }
  
  if (length(collections) == 0) {
    stop("No gene set collection selected. Please set at least one of run_hallmark/run_gobp/run_reactome to TRUE.")
  }
  
  ## -------------------------
  ## 1.4 Run GSEA
  ## 1.4 跑GSEA
  ## -------------------------
  gsea_results <- list()
  gsea_tables <- list()
  
  for (nm in names(collections)) {
    cat("Running GSEA for:", nm, "\n")
    
    term2gene <- collections[[nm]]
    
    gsea_obj <- GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      eps = eps,
      verbose = FALSE
    )
    
    gsea_results[[nm]] <- gsea_obj
    
    if (!is.null(gsea_obj) && nrow(as.data.frame(gsea_obj)) > 0) {
      gsea_tab <- as.data.frame(gsea_obj) %>%
        arrange(p.adjust)
      
      gsea_tables[[nm]] <- gsea_tab
      
      cat("  Significant/returned terms:", nrow(gsea_tab), "\n")
      
      if (save_csv) {
        out_file <- paste0(output_prefix, "_", nm, "_results.csv")
        write.csv(gsea_tab, out_file, row.names = FALSE)
        cat("  Saved:", out_file, "\n")
      }
      
    } else {
      gsea_tables[[nm]] <- data.frame()
      cat("  No enriched terms returned.\n")
    }
    
    cat("\n")
  }
  
  ## -------------------------
  ## 1.5 Return result list
  ## 1.5 返回结果
  ## -------------------------
  return(list(
    gene_list = gene_list,
    gsea_objects = gsea_results,
    gsea_tables = gsea_tables,
    collections = collections
  ))
}

#跑GSEA
gsea_out <- run_common_gsea(
  input_file = "GSE160792_DESeq2_all_results.csv",
  gene_col = "gene",
  stat_col = "stat",
  run_hallmark = TRUE,
  run_gobp = TRUE,
  run_reactome = TRUE,
  save_csv = TRUE,
  output_prefix = "GSE160792"
)


#提取geneset表
gobp_df <- gsea_out$collections$GO_BP
reactome_df <- gsea_out$collections$Reactome
hallmark_df <- gsea_out$collections$Hallmark

#提取Macrophage/NK通路
unique(gobp_df$gs_name[grepl("MACROPHAGE", gobp_df$gs_name, ignore.case = TRUE)])

unique(gobp_df$gs_name[grepl("MONOCYTE|MYELOID", gobp_df$gs_name, ignore.case = TRUE)])

unique(gobp_df$gs_name[grepl("NATURAL_KILLER|NK", gobp_df$gs_name, ignore.case = TRUE)])

unique(gobp_df$gs_name[grepl("CYTOTOXIC|CELL_KILLING", gobp_df$gs_name, ignore.case = TRUE)])

macrophage_terms <- c(
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_MACROPHAGE_DIFFERENTIATION",
  "GOBP_REGULATION_OF_MACROPHAGE_ACTIVATION"
)

nk_terms <- c(
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY"
)

macrophage_genes <- gobp_df %>%
  filter(gs_name %in% macrophage_terms) %>%
  pull(gene_symbol) %>%
  unique()

nk_genes <- gobp_df %>%
  filter(gs_name %in% nk_terms) %>%
  pull(gene_symbol) %>%
  unique()

#read DEGs file
df <- read_csv("GSE160792_DESeq2_all_results.csv", show_col_types = FALSE)

rank_df <- df %>%
  select(gene, stat) %>%
  filter(!is.na(gene), !is.na(stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(stat) %>%
  mutate(rank = row_number())

rank_df <- rank_df %>%
  mutate(
    category = case_when(
      gene %in% macrophage_genes ~ "Macrophage",
      gene %in% nk_genes ~ "NK",
      TRUE ~ "Other"
    )
  )
