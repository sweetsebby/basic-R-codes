# =========================================================
# Volcano plot script / 火山图脚本
#
# Main features / 主要功能:
# 1. Read DEG result table / 读取差异分析结果表
# 2. Symmetric x-axis / x轴左右对称
# 3. Optional automatic or manual y-axis / 可选择自动或手动设置y轴范围
# 4. Rounded background blocks / 圆角背景色块
# 5. Optional black border for Up/Down points only / 可选仅给上调和下调点加黑色边框
# 6. Manual labels, top labels, or both / 支持手动标注、前N个标注、或两者结合
# 7. Cleaner appearance for crowded points / 减轻点重叠，让图更清爽
# 8. Top count text shown above background blocks / 顶部上下调计数显示在背景框上方
# 9. P label position can be adjusted separately / P值标签位置可单独调整
# 10. Save as PNG and PDF / 保存为PNG和PDF
# =========================================================

# -------------------------
# 0. Load required packages
# 0. 加载需要的R包
# -------------------------
# install.packages(c("ggplot2", "dplyr", "ggrepel", "readr", "grid","cowplot"))

library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(grid)
library(cowplot)

# -------------------------
# 1. Input / output files
# 1. 输入与输出文件
# -------------------------
input_file  <- "GSE160792_DESeq2_all_results.csv"
output_png  <- "volcano_plot.png"
output_pdf  <- "volcano_plot.pdf"

# -------------------------
# 2. Read data first
# 2. 先读入数据
# -------------------------
df <- read_csv(input_file, show_col_types = FALSE)

# -------------------------
# 3. Column name settings
# 3. 列名设置
# Please change these names to match your DEG table exactly
# 请把下面三个列名改成与你的DEG表完全一致
# -------------------------
gene_col <- "gene"             # gene symbol column / 基因名列
fc_col   <- "log2FoldChange"   # fold change column / log2FC列
p_col    <- "pvalue"           # p value column / P值列

# -------------------------
# 4. Threshold settings
# 4. 阈值设置
# -------------------------
fc_cutoff <- 1
p_cutoff  <- 0.05

# Horizontal cutoff line on y axis
# y轴上的横向阈值线
y_cutoff_line <- -log10(p_cutoff)

# -------------------------
# 5. Axis settings
# 5. 坐标轴设置
# -------------------------
auto_axis <- TRUE         # automatic x-axis / x轴自动
auto_y_axis <- FALSE       # automatic y-axis / y轴自动

# x-axis is always symmetric around 0
# x轴始终以0为中心左右对称
x_pad <- 0.5

# Additional headroom on y-axis
# y轴额外留白
y_pad <- 1.5

# Manual axis limits if automatic mode is off
# 如果关闭自动模式，则使用手动设置
manual_x_limit <- 6
manual_y_max   <- 15

# If auto_y_axis = TRUE, whether to use quantile instead of max
# 若auto_y_axis = TRUE，是否用分位数代替最大值
use_quantile_for_auto_y <- TRUE
y_quantile <- 0.995

# -------------------------
# 6. Colors and style
# 6. 颜色与样式设置
# -------------------------
col_up   <- "#10A37F"
col_down <- "#F39C12"
col_ns   <- "#D9D9D9"

bg_left  <- "#EEDFCC"
bg_right <- "#C9D7C8"

# Point size and transparency
# 点大小和透明度
point_size  <- 2.6
point_alpha <- 0.8

# NS points smaller and lighter
# NS点更小更浅
ns_point_size  <- 1.8
ns_point_alpha <- 0.55

# Optional black border for Up/Down only
# 仅上调/下调点加黑边
add_border_updown <- TRUE
border_color <- "black"
border_size  <- 0.25

# Text sizes
# 文字大小
axis_text_size      <- 14
axis_title_size     <- 18
count_text_size     <- 7
gene_text_size      <- 4.5
p_label_text_size   <- 5
top_number_text_size <- 18
top_title_text_size  <- 18

# Output figure size
# 输出图片尺寸
out_width  <- 6
out_height <- 5
out_dpi    <- 600

# -------------------------
# 7. Label settings
# 7. 标签设置
# -------------------------

# Label mode:
# "manual" = only manual genes / 只标手动指定基因
# "top"    = only top genes / 只标前N个基因
# "both"   = manual + top / 手动和top都标
# "none"   = no labels / 不标任何基因
label_mode <- "both"

# Manually selected genes
# 手动指定要标的基因
manual_genes <- c("SPINK1", "SLC22A1", "FUT7", "CD47")

# Number of top genes to label
# 自动标注上调/下调前N个基因
top_n_up   <- 10
top_n_down <- 10

# Ranking method for top genes
# top基因排序方式
# "p"     = by smallest p value / 按最小P值
# "fc"    = by largest fold change / 按最大fold change
# "score" = by combined score / 按综合分数
top_rank_by <- "p"

# -------------------------
# 8. Check columns
# 8. 检查列名是否存在
# -------------------------
required_cols <- c(gene_col, fc_col, p_col)
missing_cols <- required_cols[!required_cols %in% colnames(df)]

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nPlease check gene_col / fc_col / p_col.\n",
      "缺少这些必要列，请检查 gene_col / fc_col / p_col 是否与表格列名一致。"
    )
  )
}

# -------------------------
# 9. Standardize columns
# 9. 统一列名并整理数据
# -------------------------
df2 <- df %>%
  transmute(
    gene   = as.character(.data[[gene_col]]),
    logFC  = as.numeric(.data[[fc_col]]),
    pvalue = as.numeric(.data[[p_col]])
  ) %>%
  filter(!is.na(gene), !is.na(logFC), !is.na(pvalue)) %>%
  mutate(
    pvalue = ifelse(pvalue <= 0, 1e-300, pvalue),
    negLog10P = -log10(pvalue),
    group = case_when(
      logFC >=  fc_cutoff & pvalue < p_cutoff ~ "Up",
      logFC <= -fc_cutoff & pvalue < p_cutoff ~ "Down",
      TRUE                                    ~ "NS"
    ),
    score = negLog10P * abs(logFC)
  )

df2$group <- factor(df2$group, levels = c("Down", "Up", "NS"))

# -------------------------
# 10. Count up/down genes
# 10. 统计上调和下调基因数
# -------------------------
n_up   <- sum(df2$group == "Up", na.rm = TRUE)
n_down <- sum(df2$group == "Down", na.rm = TRUE)

# -------------------------
# 11. Decide label genes
# 11. 决定哪些基因要标注
# -------------------------
manual_label_df <- df2 %>%
  filter(gene %in% manual_genes)

if (top_rank_by == "p") {
  top_up_df <- df2 %>%
    filter(group == "Up") %>%
    arrange(pvalue, desc(abs(logFC))) %>%
    slice_head(n = top_n_up)
  
  top_down_df <- df2 %>%
    filter(group == "Down") %>%
    arrange(pvalue, desc(abs(logFC))) %>%
    slice_head(n = top_n_down)
  
} else if (top_rank_by == "fc") {
  top_up_df <- df2 %>%
    filter(group == "Up") %>%
    arrange(desc(logFC), pvalue) %>%
    slice_head(n = top_n_up)
  
  top_down_df <- df2 %>%
    filter(group == "Down") %>%
    arrange(logFC, pvalue) %>%
    slice_head(n = top_n_down)
  
} else if (top_rank_by == "score") {
  top_up_df <- df2 %>%
    filter(group == "Up") %>%
    arrange(desc(score), pvalue) %>%
    slice_head(n = top_n_up)
  
  top_down_df <- df2 %>%
    filter(group == "Down") %>%
    arrange(desc(score), pvalue) %>%
    slice_head(n = top_n_down)
  
} else {
  stop("top_rank_by must be one of: 'p', 'fc', 'score'")
}

top_label_df <- bind_rows(top_up_df, top_down_df)

if (label_mode == "manual") {
  label_df <- manual_label_df
} else if (label_mode == "top") {
  label_df <- top_label_df
} else if (label_mode == "both") {
  label_df <- bind_rows(manual_label_df, top_label_df) %>%
    distinct(gene, .keep_all = TRUE)
} else if (label_mode == "none") {
  label_df <- df2[0, ]
} else {
  stop("label_mode must be one of: 'manual', 'top', 'both', 'none'")
}

# -------------------------
# 12. Axis range
# 12. 坐标轴范围设置
# -------------------------

# ---- x-axis / x轴 ----
if (auto_axis) {
  x_abs_max <- max(abs(df2$logFC), na.rm = TRUE)
  x_limit <- ceiling(x_abs_max + x_pad)
  x_min <- -x_limit
  x_max <-  x_limit
} else {
  x_min <- -manual_x_limit
  x_max <-  manual_x_limit
}

# ---- y-axis / y轴 ----
if (auto_y_axis) {
  if (use_quantile_for_auto_y) {
    y_main_max <- quantile(df2$negLog10P, y_quantile, na.rm = TRUE)
  } else {
    y_main_max <- max(df2$negLog10P, na.rm = TRUE)
  }
  
  y_label_max <- if (nrow(label_df) > 0) max(label_df$negLog10P, na.rm = TRUE) else y_main_max
  
  # 数据区上限
  # Upper limit for data region only
  y_max <- ceiling(max(y_main_max, min(y_label_max, y_main_max + 2)) + y_pad)
  
} else {
  y_max <- manual_y_max
}

# -------------------------
# 13. Background settings
# 13. 背景色块设置
# -------------------------
# 背景框只跟数据区高度走
# Background block follows data-region height only
bg_ymin <- y_cutoff_line
bg_ymax <- y_max * 0.96

left_xmin  <- x_min + 0.5
left_xmax  <- -0.2
right_xmin <- 0.2
right_xmax <- x_max - 0.5

x_to_npc <- function(x, xmin, xmax) {
  (x - xmin) / (xmax - xmin)
}

y_to_npc <- function(y, ymin, ymax) {
  (y - ymin) / (ymax - ymin)
}

left_center_x  <- (left_xmin + left_xmax) / 2
right_center_x <- (right_xmin + right_xmax) / 2
block_center_y <- (bg_ymin + bg_ymax) / 2

left_width   <- (left_xmax - left_xmin) / (x_max - x_min)
right_width  <- (right_xmax - right_xmin) / (x_max - x_min)
block_height <- (bg_ymax - bg_ymin) / (y_max - 0)

# -------------------------
# 13.5 Text position settings
# 13.5 文字位置设置
# -------------------------

# Only x positions are set manually
# 这里只手动设置左右横坐标
top_text_x_left  <- x_min * 0.48
top_text_x_right <- x_max * 0.48

# Top text positions are computed from current y_max,
# but they do NOT change the background height
# 顶部数字/标题位置根据当前y_max自动计算，
# 但不会反过来影响背景框高度
top_title_y  <- y_max + max(1.2, 0.06 * y_max)
top_number_y <- y_max + max(3.0, 0.14 * y_max)

# Final visible top limit
# 最终整张图的显示上限（包含顶部文字）
plot_top_y <- top_number_y + max(1, 0.03 * y_max)

# P label position
# P值标签位置
p_label_offset <- max(0.35, 0.03 * y_max)
p_label_x <- x_max * 0.78
p_label_y <- y_cutoff_line + p_label_offset

# -------------------------
# 14. Plot
# 14. 作图
# -------------------------
p <- ggplot(df2, aes(x = logFC, y = negLog10P)) +
  
  # Rounded background blocks
  # 圆角背景色块
  annotation_custom(
    grob = roundrectGrob(
      x = unit(x_to_npc(left_center_x, x_min, x_max), "npc"),
      y = unit(y_to_npc(block_center_y, 0, y_max), "npc"),
      width = unit(left_width, "npc"),
      height = unit(block_height, "npc"),
      r = unit(0.035, "snpc"),
      gp = gpar(fill = bg_left, col = NA, alpha = 0.45)
    )
  ) +
  annotation_custom(
    grob = roundrectGrob(
      x = unit(x_to_npc(right_center_x, x_min, x_max), "npc"),
      y = unit(y_to_npc(block_center_y, 0, y_max), "npc"),
      width = unit(right_width, "npc"),
      height = unit(block_height, "npc"),
      r = unit(0.035, "snpc"),
      gp = gpar(fill = bg_right, col = NA, alpha = 0.45)
    )
  ) +
  
  # NS points
  # NS点
  geom_point(
    data = df2 %>% filter(group == "NS"),
    color = col_ns,
    size = ns_point_size,
    alpha = ns_point_alpha,
    stroke = 0
  ) +
  
  # Down-regulated points
  # 下调点
  geom_point(
    data = df2 %>% filter(group == "Down"),
    shape = 21,
    fill = col_down,
    color = if (add_border_updown) border_color else col_down,
    stroke = if (add_border_updown) border_size else 0,
    size = point_size,
    alpha = point_alpha
  ) +
  
  # Up-regulated points
  # 上调点
  geom_point(
    data = df2 %>% filter(group == "Up"),
    shape = 21,
    fill = col_up,
    color = if (add_border_updown) border_color else col_up,
    stroke = if (add_border_updown) border_size else 0,
    size = point_size,
    alpha = point_alpha
  ) +
  
  # Cutoff lines
  # 阈值线
  geom_hline(
    yintercept = y_cutoff_line,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # Gene labels
  # 基因标签
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = gene_text_size,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = NA,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  
  
  # P-value label
  # P值标签
  annotate(
    "text",
    x = p_label_x,
    y = p_label_y,
    label = paste0("P = ", p_cutoff),
    size = p_label_text_size,
    fontface = "italic"
  ) +
  
  # Axis scales
  # 坐标轴刻度
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = pretty(c(x_min, x_max), n = 7)
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = pretty(c(0, y_max), n = 5),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  # Axis labels
  # 坐标轴标题
  labs(
    x = "log2(FC)",
    y = expression(-log[10](italic(P)~value))
  ) +
  
  # Theme
  # 主题样式
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_title_size, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    plot.margin = margin(0, 25, 20, 25),
    panel.background = element_rect(fill = "white", color = NA)
  )



# -------------------------
# 15. Add top labels outside plot panel
# 15. 在图外上方添加顶部标注
# -------------------------

top_layer <- ggdraw() +
  draw_label(
    label = as.character(n_down),
    x = 0.35, y = 0.58,
    vjust = 0.5,
    size = top_number_text_size
  ) +
  draw_label(
    label = "Down-regulated",
    x = 0.35, y = 0.02,
    vjust = -0.3,
    size = top_title_text_size
  ) +
  draw_label(
    label = as.character(n_up),
    x = 0.70, y = 0.58,
    vjust = 0.5,
    size = top_number_text_size
  ) +
  draw_label(
    label = "Up-regulated",
    x = 0.70, y = 0.02,
    vjust = -0.3,
    size = top_title_text_size
  )

final_plot <- plot_grid(
  top_layer,
  p,
  ncol = 1,
  rel_heights = c(0.12, 1)
)
# -------------------------
# 15. Show plot
# 15. 显示图像
# -------------------------
print(final_plot)

# -------------------------
# 16. Save files
# 16. 保存图片
# -------------------------
ggsave(
  filename = output_png,
  plot = final_plot,
  width = out_width,
  height = out_height,
  dpi = out_dpi,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = out_width,
  height = out_height,
  bg = "white"
)