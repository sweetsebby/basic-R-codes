# =========================================================
# Rank plot - manual only version
# 纯手动版 rank plot
# =========================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(grid)

# =========================================================
# 0. Read DEG result and build rank_df
# 0. 读入DEG结果并生成 rank_df
# =========================================================
df <- read_csv("GSE160792_DESeq2_all_results.csv", show_col_types = FALSE)

rank_df <- df %>%
  select(gene, stat) %>%
  filter(!is.na(gene), !is.na(stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(stat) %>%
  mutate(rank = row_number())

# =========================================================
# 1. Manual gene groups
# 1. 手动指定基因分组
# =========================================================
# 这里直接手动指定哪些属于 macrophage，哪些属于 NK
# Directly define which genes belong to macrophage or NK

manual_macrophage_genes <- c("CD47", "ICAM1", "CD84", "SYK", "BPI")
manual_nk_genes <- c("NKG7", "PRF1", "GZMB", "NFIL3", "RAB27A")

# 所有手动标签基因
# All manually labeled genes
manual_genes <- unique(c(manual_macrophage_genes, manual_nk_genes))

# =========================================================
# 2. Add category
# 2. 添加类别信息
# =========================================================
rank_df <- rank_df %>%
  mutate(
    category = case_when(
      gene %in% manual_macrophage_genes ~ "Macrophage",
      gene %in% manual_nk_genes ~ "NK",
      TRUE ~ "Other"
    )
  )

plot_df <- rank_df

label_df <- plot_df %>%
  filter(gene %in% manual_genes)

# =========================================================
# 3. Plot appearance settings
# 3. 作图外观参数
# =========================================================

# colors / 颜色
col_other <- "grey70"
col_mac   <- "#E69F00"
col_nk    <- "#8E0D8A"

fill_mac  <- "#F3D68A"
fill_nk   <- "#D9A3D9"

# point sizes / 点大小
other_point_size  <- 1.0
group_point_size  <- 2.2
legend_point_size <- 4

# label sizes / 标签大小
label_text_size   <- 4.2
segment_line_size <- 0.45

# axis / legend text / 坐标轴和图例文字大小
axis_text_size    <- 18
axis_title_size   <- 24
legend_text_size  <- 15

# axis labels / 坐标轴标题
x_label <- "Rank"
y_label <- "NormZ"   # 更严谨可以改成 "stat"

# top-left annotation / 左上角文字
top_left_label <- "5-aza-dC"
top_left_x <- 0
top_left_y <- 25
top_left_text_size <- 8

# legend position / 图例位置
legend_x <- 0.80
legend_y <- 0.10

# output settings / 输出设置
output_png <- "rank_plot_manual.png"
output_pdf <- "rank_plot_manual.pdf"
out_width  <- 8
out_height <- 8
out_dpi    <- 600

# =========================================================
# 4. Plot
# 4. 作图
# =========================================================
p_rank <- ggplot(plot_df, aes(x = rank, y = stat)) +
  
  # all genes in grey / 所有基因灰点
  geom_point(
    data = plot_df,
    color = col_other,
    size = other_point_size,
    alpha = 1
  ) +
  
  # manually selected macrophage genes / 手动指定的 macrophage 基因
  geom_point(
    data = label_df %>% filter(category == "Macrophage"),
    color = col_mac,
    size = group_point_size
  ) +
  
  # manually selected NK genes / 手动指定的 NK 基因
  geom_point(
    data = label_df %>% filter(category == "NK"),
    color = col_nk,
    size = group_point_size
  ) +
  
  # labels with colored boxes / 带底色标签框
  geom_label_repel(
    data = label_df,
    aes(label = gene, fill = category),
    color = "black",
    label.size = NA,
    size = label_text_size,
    box.padding = 0.20,
    point.padding = 0.15,
    segment.color = "black",
    segment.size = segment_line_size,
    segment.alpha = 1,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  
  # label fill colors / 标签底色
  scale_fill_manual(
    values = c(
      "Macrophage" = fill_mac,
      "NK" = fill_nk
    )
  ) +
  
  # legend colors / 图例颜色
  scale_color_manual(
    values = c(
      "Macrophage" = col_mac,
      "NK" = col_nk
    )
  ) +
  
  # invisible points for legend / 用虚拟点生成图例
  geom_point(
    data = data.frame(
      rank = c(NA, NA),
      stat = c(NA, NA),
      category = c("Macrophage", "NK")
    ),
    aes(color = category),
    size = legend_point_size
  ) +
  
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = legend_point_size + 1),
      nrow = 2,
      byrow = TRUE
    )
  ) +
  
  labs(
    x = x_label,
    y = y_label
  ) +
  
  # top-left text / 左上角文字
  annotate(
    "text",
    x = top_left_x,
    y = top_left_y,
    label = top_left_label,
    hjust = 0,
    vjust = 1,
    size = top_left_text_size
  ) +
  
  theme_bw(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_title_size, color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = unit(0.22, "cm"),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    legend.position = c(legend_x, legend_y),
    legend.direction = "vertical",
    legend.text = element_text(size = legend_text_size, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    
    plot.margin = margin(20, 20, 20, 20)
  )

# show plot / 显示图
print(p_rank)

# =========================================================
# 5. Save plot
# 5. 保存图片
# =========================================================
ggsave(
  filename = output_png,
  plot = p_rank,
  width = out_width,
  height = out_height,
  dpi = out_dpi,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p_rank,
  width = out_width,
  height = out_height,
  bg = "white"
)