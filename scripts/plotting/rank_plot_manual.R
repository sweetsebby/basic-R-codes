# =========================================================
# Publication-style rank plot (2 groups, manual labels)
# 发表风格版 rank plot（双组，手动标签）
# =========================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(grid)

# =========================================================
# 0. Input and DEG table settings
# 0. 输入文件与列名设置
# =========================================================
input_file <- "GSE160792_DESeq2_all_results.csv"

gene_col <- "gene"
stat_col <- "stat"

# =========================================================
# 1. Read data and build rank_df
# 1. 读取数据并生成 rank_df
# =========================================================
df <- read_csv(input_file, show_col_types = FALSE)

rank_df <- df %>%
  transmute(
    gene = as.character(.data[[gene_col]]),
    stat = as.numeric(.data[[stat_col]])
  ) %>%
  filter(!is.na(gene), !is.na(stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(stat) %>%
  mutate(rank = row_number())

# =========================================================
# 2. Manual two-group definition
# 2. 手动双组定义
# =========================================================
group1_name <- "Group1"
group2_name <- "Group2"

group1_genes <- c("CD47", "ICAM1", "CD84", "SYK", "BPI")
group2_genes <- c("NKG7", "PRF1", "GZMB", "NFIL3", "RAB27A")

# =========================================================
# 3. Plot labels and text settings
# 3. 图中文字设置
# =========================================================
top_left_label <- "example"
top_left_x <- 0
top_left_y <- 25
top_left_text_size <- 9

x_label <- "Rank"
y_label <- "DESeq2 statistic"

# =========================================================
# 4. Style settings
# 4. 样式设置
# =========================================================

# ---- colors / 颜色 ----
col_other  <- "#D9D9D9"
col_group1 <- "#10A37F"
col_group2 <- "#F39C12"

fill_group1 <- "#C3D7C2"
fill_group2 <- "#EBCB8F"

# ---- point sizes / 点大小 ----
other_point_size  <- 0.9
group_point_size  <- 2.3
legend_point_size <- 4.2

# ---- label sizes / 标签大小 ----
label_text_size     <- 4.0
segment_line_size   <- 0.32
label_box_padding   <- 0.22
label_point_padding <- 0.18

# ---- axis / legend / 坐标轴与图例 ----
axis_text_size   <- 18
axis_title_size  <- 24
legend_text_size <- 15

# ---- legend position / 图例位置 ----
legend_x <- 0.81
legend_y <- 0.12

# ---- label behavior / 标签排斥参数 ----
repel_force <- 1.8
repel_force_pull <- 0.35

# ---- output / 输出 ----
output_png <- "rank_plot_pubstyle_2groups.png"
output_pdf <- "rank_plot_pubstyle_2groups.pdf"
out_width  <- 14
out_height <- 8
out_dpi    <- 600

# =========================================================
# 5. Add category information
# 5. 添加类别信息
# =========================================================
plot_df <- rank_df %>%
  mutate(
    category = case_when(
      gene %in% group1_genes ~ group1_name,
      gene %in% group2_genes ~ group2_name,
      TRUE ~ "Other"
    )
  )

label_df <- plot_df %>%
  filter(category != "Other")

plot_df$category <- factor(plot_df$category, levels = c(group1_name, group2_name, "Other"))
label_df$category <- factor(label_df$category, levels = c(group1_name, group2_name))

# =========================================================
# 6. Axis helper values
# 6. 坐标辅助参数
# =========================================================
x_max <- max(plot_df$rank, na.rm = TRUE)
y_min <- min(plot_df$stat, na.rm = TRUE)
y_max <- max(plot_df$stat, na.rm = TRUE)

top_padding <- 0.05 * (y_max - y_min)

# color maps / 颜色映射
fill_map <- c(fill_group1, fill_group2)
names(fill_map) <- c(group1_name, group2_name)

color_map <- c(col_group1, col_group2)
names(color_map) <- c(group1_name, group2_name)

# legend helper dataframe / 图例辅助数据
legend_df <- data.frame(
  rank = c(NA, NA),
  stat = c(NA, NA),
  category = factor(c(group1_name, group2_name), levels = c(group1_name, group2_name))
)

# =========================================================
# 7. Plot
# 7. 作图
# =========================================================
p_rank <- ggplot(plot_df, aes(x = rank, y = stat)) +
  
  geom_point(
    color = col_other,
    size = other_point_size,
    alpha = 0.9
  ) +
  
  geom_point(
    data = label_df %>% filter(category == group1_name),
    color = col_group1,
    size = group_point_size
  ) +
  
  geom_point(
    data = label_df %>% filter(category == group2_name),
    color = col_group2,
    size = group_point_size
  ) +
  
  geom_label_repel(
    data = label_df,
    aes(label = gene, fill = category),
    color = "black",
    label.size = 0.2,
    label.r = unit(0.12, "lines"),
    size = label_text_size,
    box.padding = label_box_padding,
    point.padding = label_point_padding,
    segment.color = "grey25",
    segment.size = segment_line_size,
    segment.alpha = 0.9,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = repel_force,
    force_pull = repel_force_pull,
    seed = 123,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = fill_map) +
  
  geom_point(
    data = legend_df,
    aes(color = category),
    size = legend_point_size
  ) +
  
  scale_color_manual(values = color_map) +
  
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = legend_point_size + 0.8),
      ncol = 1,
      byrow = TRUE
    )
  ) +
  
  labs(
    x = x_label,
    y = y_label
  ) +
  
  annotate(
    "text",
    x = top_left_x,
    y = top_left_y,
    label = top_left_label,
    hjust = 0,
    vjust = 1,
    size = top_left_text_size
  ) +
  
  coord_cartesian(
    xlim = c(0, x_max * 1.02),
    ylim = c(y_min, y_max + top_padding),
    clip = "on"
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_title_size, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.9),
    axis.ticks = element_line(color = "black", linewidth = 0.9),
    axis.ticks.length = unit(0.20, "cm"),
    legend.position = c(legend_x, legend_y),
    legend.direction = "vertical",
    legend.text = element_text(size = legend_text_size, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(18, 22, 18, 18)
  )

# =========================================================
# 8. Show plot
# 8. 显示图
# =========================================================
print(p_rank)

# =========================================================
# 9. Save plot
# 9. 保存图片
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
