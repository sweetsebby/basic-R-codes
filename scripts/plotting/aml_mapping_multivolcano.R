# =========================================================
# GSE116256 AML cluster-to-healthy mapping and core-cluster volcano plots
# 完整整理版脚本（统一配色 + 拼图输出）
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scRNAtoolVis)
  library(corrplot)
  library(grid)
})

# =========================================================
# 0. Global settings
# =========================================================

base_out_dir <- "~/AML_scRNA_project"
out_dir_map  <- file.path(base_out_dir, "healthy_reference_mapping")
out_dir_mv   <- file.path(base_out_dir, "core_cluster_multivolcano")
out_dir_comb <- file.path(base_out_dir, "combined_figures")

dir.create(out_dir_map,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_mv,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_comb, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 0.1 Unified color palette (your RGB colors)
# 统一使用你给的RGB配色
# =========================================================

my_cols8 <- c(
  rgb(29,  82, 161, maxColorValue = 255),
  rgb(113,109, 178, maxColorValue = 255),
  rgb(101,200, 204, maxColorValue = 255),
  rgb(114,193,  90, maxColorValue = 255),
  rgb(240,233,  75, maxColorValue = 255),
  rgb(243,121,  59, maxColorValue = 255),
  rgb(198, 81, 159, maxColorValue = 255),
  rgb(243,135, 141, maxColorValue = 255)
)

# =========================================================
# 1. Load processed objects
# =========================================================

aml_union   <- readRDS("~/AML_scRNA_project/GSE116256_core_signature_full_analysis/OBJECT_aml_union_processed.rds")
healthy_obj <- readRDS("~/AML_scRNA_project/GSE116256_core_signature_full_analysis/AML_niche_interaction_analysis/GSE116256_Healthy_rebuilt_processed.rds")

# Quick checks
print(Reductions(aml_union))
print(Reductions(healthy_obj))

print(table(aml_union$seurat_clusters, useNA = "ifany"))
print(table(healthy_obj$seurat_clusters, useNA = "ifany"))

# =========================================================
# 2. Re-annotate healthy reference
# =========================================================

healthy_obj$major_celltype <- as.character(healthy_obj$seurat_clusters)

healthy_obj$major_celltype[healthy_obj$seurat_clusters == "0"]  <- "HSPC_lymphoid_primed"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "1"]  <- "HSPC_stem_like"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "2"]  <- "Monocyte_macrophage_like"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "3"]  <- "T_cell"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "4"]  <- "Erythroid"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "5"]  <- "Inflammatory_monocyte"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "6"]  <- "Neutrophil_granulocytic_precursor"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "7"]  <- "Pre_B_cell"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "8"]  <- "NK_cytotoxic"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "9"]  <- "B_cell_mature"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "10"] <- "Plasma_cell"
healthy_obj$major_celltype[healthy_obj$seurat_clusters == "11"] <- "pDC"

healthy_obj$major_celltype <- factor(healthy_obj$major_celltype)

healthy_obj$broad_cellclass <- NA_character_

healthy_obj$broad_cellclass[healthy_obj$major_celltype %in% c(
  "HSPC_lymphoid_primed", "HSPC_stem_like"
)] <- "HSPC"

healthy_obj$broad_cellclass[healthy_obj$major_celltype %in% c(
  "Monocyte_macrophage_like", "Inflammatory_monocyte", "Neutrophil_granulocytic_precursor", "pDC"
)] <- "Myeloid"

healthy_obj$broad_cellclass[healthy_obj$major_celltype %in% c(
  "T_cell", "NK_cytotoxic"
)] <- "Lymphoid_cytotoxic"

healthy_obj$broad_cellclass[healthy_obj$major_celltype %in% c(
  "Pre_B_cell", "B_cell_mature", "Plasma_cell"
)] <- "B_lineage"

healthy_obj$broad_cellclass[healthy_obj$major_celltype %in% c(
  "Erythroid"
)] <- "Erythroid"

healthy_obj$broad_cellclass <- factor(healthy_obj$broad_cellclass)

print(table(healthy_obj$major_celltype, useNA = "ifany"))
print(table(healthy_obj$broad_cellclass, useNA = "ifany"))

saveRDS(healthy_obj, file.path(out_dir_map, "healthy_obj_reannotated.rds"))

# =========================================================
# 3. Label transfer: healthy reference -> AML query
# =========================================================

DefaultAssay(aml_union)   <- "RNA"
DefaultAssay(healthy_obj) <- "RNA"

shared_genes <- intersect(rownames(aml_union), rownames(healthy_obj))
cat("Shared genes:", length(shared_genes), "\n")

aml_query   <- subset(aml_union, features = shared_genes)
healthy_ref <- subset(healthy_obj, features = shared_genes)

anchors <- FindTransferAnchors(
  reference = healthy_ref,
  query = aml_query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)

pred_major <- TransferData(
  anchorset = anchors,
  refdata = healthy_ref$major_celltype,
  dims = 1:30
)

aml_query <- AddMetaData(aml_query, metadata = pred_major)

print(table(aml_query$predicted.id, useNA = "ifany"))
print(summary(aml_query$prediction.score.max))

saveRDS(aml_query, file.path(out_dir_map, "aml_query_with_transferred_labels.rds"))

# =========================================================
# 4. Prepare color mapping for predicted.id
# =========================================================

pred_levels <- levels(factor(aml_query$predicted.id))
pred_cols   <- colorRampPalette(my_cols8)(length(pred_levels))
names(pred_cols) <- pred_levels

# =========================================================
# 5. UMAP visualization
# =========================================================

# 5.1 UMAP colored by predicted identity
p_pred <- DimPlot(
  aml_query,
  reduction = "umap",
  group.by = "predicted.id",
  label = TRUE,
  repel = TRUE,
  cols = pred_cols
) +
  ggtitle("Transferred healthy-like identity") +
  theme_bw(base_size = 12) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave(
  file.path(out_dir_map, "AML_query_UMAP_predicted_identity.png"),
  p_pred, width = 10, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "AML_query_UMAP_predicted_identity.pdf"),
  p_pred, width = 10, height = 8, bg = "white"
)

# 5.2 UMAP colored by prediction score
p_score <- FeaturePlot(
  aml_query,
  features = "prediction.score.max",
  reduction = "umap"
) +
  scale_colour_gradientn(colours = c("white", my_cols8[3], my_cols8[1])) +
  ggtitle("Label transfer prediction score") +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir_map, "AML_query_UMAP_prediction_score.png"),
  p_score, width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "AML_query_UMAP_prediction_score.pdf"),
  p_score, width = 8, height = 6, bg = "white"
)

# 5.3 Cluster vs predicted identity side-by-side
p_cluster <- DimPlot(
  aml_query,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  scale_color_manual(values = colorRampPalette(my_cols8)(length(unique(aml_query$seurat_clusters)))) +
  ggtitle("AML clusters") +
  theme_bw(base_size = 12)

p_cluster_pred_side <- p_cluster + p_pred

ggsave(
  file.path(out_dir_map, "AML_query_UMAP_cluster_vs_prediction_side_by_side.png"),
  p_cluster_pred_side, width = 16, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "AML_query_UMAP_cluster_vs_prediction_side_by_side.pdf"),
  p_cluster_pred_side, width = 16, height = 8, bg = "white"
)

# =========================================================
# 6. Cluster-level summary tables
# =========================================================

cluster_pred_tbl <- aml_query@meta.data %>%
  as.data.frame() %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(prediction.score.max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(
  cluster_pred_tbl,
  file.path(out_dir_map, "AML_cluster_by_predicted_identity.csv")
)

cluster_top_tbl <- cluster_pred_tbl %>%
  group_by(seurat_clusters) %>%
  arrange(desc(frac), desc(mean_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write_csv(
  cluster_top_tbl,
  file.path(out_dir_map, "AML_cluster_top_predicted_identity.csv")
)

print(cluster_top_tbl)

# =========================================================
# 7. Heatmap of cluster × predicted identity
# =========================================================

p_heat <- ggplot(cluster_pred_tbl, aes(x = predicted.id, y = factor(seurat_clusters), fill = frac)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("white", my_cols8[4], my_cols8[6])) +
  theme_bw(base_size = 12) +
  labs(
    x = "Predicted healthy BM cell type",
    y = "AML cluster",
    fill = "Fraction"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir_map, "AML_cluster_by_predicted_identity_heatmap.png"),
  p_heat, width = 10, height = 6, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "AML_cluster_by_predicted_identity_heatmap.pdf"),
  p_heat, width = 10, height = 6, bg = "white"
)

# =========================================================
# 8. Focus on core clusters (2/5/9/10)
# =========================================================

core_clusters <- c("2", "5", "9", "10")

core_cluster_pred_tbl <- cluster_pred_tbl %>%
  filter(seurat_clusters %in% core_clusters) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = core_clusters))

write_csv(
  core_cluster_pred_tbl,
  file.path(out_dir_map, "core_clusters_predicted_identity_summary.csv")
)

p_core_stack <- ggplot(core_cluster_pred_tbl, aes(x = seurat_clusters, y = frac, fill = predicted.id)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pred_cols) +
  theme_bw(base_size = 12) +
  labs(
    title = "Healthy BM counterpart mapping of core-enriched AML clusters",
    x = "AML cluster",
    y = "Fraction",
    fill = "Predicted healthy-like identity"
  )

ggsave(
  file.path(out_dir_map, "core_clusters_predicted_identity_stackedbar.png"),
  p_core_stack, width = 8, height = 5.5, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "core_clusters_predicted_identity_stackedbar.pdf"),
  p_core_stack, width = 8, height = 5.5, bg = "white"
)

# =========================================================
# 9. Mark non-malignant-like clusters
# =========================================================

aml_query$cell_compartment <- "AML_like"
aml_query$cell_compartment[aml_query$seurat_clusters %in% c("0","7","8","11","14","15","17")] <- "non_malignant_like"

print(table(aml_query$cell_compartment, aml_query$predicted.id))

write_csv(
  aml_query@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id"),
  file.path(out_dir_map, "aml_query_metadata_with_prediction_and_compartment.csv")
)

cluster_top_tbl_aml_like <- aml_query@meta.data %>%
  as.data.frame() %>%
  filter(cell_compartment == "AML_like") %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(prediction.score.max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  arrange(desc(frac), desc(mean_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write_csv(
  cluster_top_tbl_aml_like,
  file.path(out_dir_map, "AML_like_cluster_top_predicted_identity.csv")
)

# =========================================================
# 10. UMAP with core-cluster labels only
# =========================================================

short_label_map <- c(
  "Inflammatory_monocyte" = "Inflamm_mono",
  "Monocyte_macrophage_like" = "Mono_macro",
  "HSPC_lymphoid_primed" = "HSPC_like",
  "HSPC_stem_like" = "HSPC_stem",
  "Neutrophil_granulocytic_precursor" = "Neut_gran_prec",
  "NK_cytotoxic" = "NK",
  "T_cell" = "T",
  "B_cell_mature" = "B",
  "Plasma_cell" = "Plasma",
  "Erythroid" = "Ery",
  "pDC" = "pDC",
  "Pre_B_cell" = "Pre_B"
)

cluster_top_tbl2 <- aml_query@meta.data %>%
  as.data.frame() %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  arrange(desc(frac), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    predicted_short = recode(predicted.id, !!!short_label_map),
    label = paste0(seurat_clusters, " | ", predicted_short)
  )

emb <- Embeddings(aml_query, "umap") %>% as.data.frame()
emb$cell <- rownames(emb)

meta_plot <- aml_query@meta.data %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))

plot_df <- left_join(emb, meta_plot, by = "cell")

cluster_centers <- plot_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    x = median(umap_1),
    y = median(umap_2),
    .groups = "drop"
  ) %>%
  left_join(cluster_top_tbl2, by = "seurat_clusters")

cluster_centers_core <- cluster_centers %>%
  filter(seurat_clusters %in% c("2", "5", "9", "10")) %>%
  mutate(
    x_label = case_when(
      seurat_clusters == "2"  ~ x + 0.2,
      seurat_clusters == "5"  ~ x + 0.2,
      seurat_clusters == "9"  ~ x - 0.8,
      seurat_clusters == "10" ~ x + 0.9,
      TRUE ~ x
    ),
    y_label = case_when(
      seurat_clusters == "2"  ~ y - 2.0,
      seurat_clusters == "5"  ~ y + 2.7,
      seurat_clusters == "9"  ~ y + 2.9,
      seurat_clusters == "10" ~ y + 2.1,
      TRUE ~ y
    )
  )

p_core_label <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = predicted.id)) +
  geom_point(size = 0.25, alpha = 0.85) +
  geom_segment(
    data = cluster_centers_core,
    aes(x = x, y = y, xend = x_label, yend = y_label),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.4
  ) +
  geom_label(
    data = cluster_centers_core,
    aes(x = x_label, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  scale_color_manual(values = pred_cols) +
  theme_bw(base_size = 12) +
  labs(
    title = "AML UMAP with top predicted identity for core clusters",
    color = "Predicted healthy-like identity"
  ) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.8, "cm")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4, alpha = 1))
  )

ggsave(
  file.path(out_dir_map, "AML_query_UMAP_core_cluster_labels.png"),
  p_core_label, width = 14, height = 9, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_map, "AML_query_UMAP_core_cluster_labels.pdf"),
  p_core_label, width = 14, height = 9, bg = "white"
)

# =========================================================
# 11. Multi-volcano plot for core AML clusters (2/5/9/10)
# =========================================================

aml_like_clusters <- c("1","2","3","4","5","6","9","10","12","13","16","18","19")
core_clusters     <- c("2","5","9","10")

aml_mv <- subset(aml_union, subset = seurat_clusters %in% aml_like_clusters)
Idents(aml_mv) <- "seurat_clusters"

deg_list_mv <- lapply(core_clusters, function(cl) {
  message("Running DEG for cluster ", cl, " vs other AML-like clusters ...")
  
  res <- FindMarkers(
    object = aml_mv,
    ident.1 = cl,
    ident.2 = setdiff(aml_like_clusters, cl),
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0
  )
  
  res %>%
    tibble::rownames_to_column("gene") %>%
    mutate(cluster = paste0("Cluster_", cl))
})

names(deg_list_mv) <- core_clusters
deg_mv_tbl <- bind_rows(deg_list_mv)

write_csv(
  deg_mv_tbl,
  file.path(out_dir_mv, "core_clusters_2_5_9_10_vs_other_AMLlike_DEG.csv")
)

# =========================================================
# 12. jjVolcano topGene mode
# =========================================================

cluster_order_raw <- c("Cluster_2", "Cluster_5", "Cluster_9", "Cluster_10")

cluster_label_map <- c(
  "Cluster_2"  = "C2 | Inflam-mono",
  "Cluster_5"  = "C5 | Mono-macro",
  "Cluster_9"  = "C9 | Mono-macro",
  "Cluster_10" = "C10 | Mono/Inflam"
)

diffData_plot <- deg_mv_tbl %>%
  filter(cluster %in% cluster_order_raw) %>%
  mutate(cluster = recode(cluster, !!!cluster_label_map))

cluster_order_plot <- unname(cluster_label_map[cluster_order_raw])

p_mv_top <- jjVolcano(
  diffData = diffData_plot,
  topGeneN = 6,
  log2FC.cutoff = 0.5,
  col.type = "updown",
  aesCol = c(my_cols8[1], my_cols8[6]),
  tile.col = my_cols8[c(1,2,3,4)],
  cluster.order = cluster_order_plot,
  size = 3.8,
  fontface = "italic",
  legend.position = c(0.96, 0.985)
)

ggsave(
  file.path(out_dir_mv, "core_clusters_multivolcano_jjVolcano_topGene_clean.png"),
  p_mv_top, width = 15, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_mv, "core_clusters_multivolcano_jjVolcano_topGene_clean.pdf"),
  p_mv_top, width = 15, height = 8, bg = "white"
)

# =========================================================
# 13. Combine UMAP + volcano into one figure
#    不加字母，直接拼一起
# =========================================================

p_combined <- p_core_label + p_mv_top +
  plot_layout(widths = c(1.0, 1.35))

ggsave(
  file.path(out_dir_comb, "AML_UMAP_and_multivolcano_combined.png"),
  p_combined, width = 22, height = 9, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir_comb, "AML_UMAP_and_multivolcano_combined.pdf"),
  p_combined, width = 22, height = 9, bg = "white"
)

# =========================================================
# 14. Save key R objects
# =========================================================

saveRDS(healthy_obj,              file.path(out_dir_map, "OBJECT_healthy_obj_reannotated.rds"))
saveRDS(healthy_ref,              file.path(out_dir_map, "OBJECT_healthy_ref_sharedgenes.rds"))
saveRDS(aml_query,                file.path(out_dir_map, "OBJECT_aml_query_with_mapping.rds"))
saveRDS(cluster_pred_tbl,         file.path(out_dir_map, "OBJECT_cluster_pred_tbl.rds"))
saveRDS(cluster_top_tbl,          file.path(out_dir_map, "OBJECT_cluster_top_tbl.rds"))
saveRDS(cluster_top_tbl_aml_like, file.path(out_dir_map, "OBJECT_cluster_top_tbl_aml_like.rds"))
saveRDS(deg_mv_tbl,               file.path(out_dir_mv,  "OBJECT_deg_mv_tbl.rds"))

message("All analysis finished successfully.")
message("Mapping results saved to: ", normalizePath(out_dir_map))
message("Volcano results saved to: ", normalizePath(out_dir_mv))
message("Combined figure saved to: ", normalizePath(out_dir_comb))
