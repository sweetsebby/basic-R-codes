# =========================================================
# Bubble plot from simple X-Y-value CSV
# 用简单 X-Y-value CSV 画气泡图
# 支持：
# 1. 多行矩阵模式
# 2. 单行模式
# 3. 自动判断模式
# 4. 蓝色马卡龙大跨度渐变
# =========================================================

library(ggplot2)
library(dplyr)
library(readr)
library(grid)

# =========================================================
# 1. Read simple bubble data
# 1. 读取简单格式数据
# =========================================================
# CSV 必须包含三列：
# X, Y, Relative protein level
#
# 例如：
# X,Y,Relative protein level
# 1,1,1
# 2,1,1
# 3,1,1
# ...

read_xy_bubble_data <- function(input_file) {
  raw_df <- read_csv(input_file, show_col_types = FALSE)
  
  required_cols <- c("X", "Y", "Relative protein level")
  missing_cols <- required_cols[!required_cols %in% colnames(raw_df)]
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing columns: ",
        paste(missing_cols, collapse = ", "),
        "\nCSV 必须包含这三列: X, Y, Relative protein level"
      )
    )
  }
  
  df <- raw_df %>%
    transmute(
      x = suppressWarnings(as.numeric(X)),
      y = suppressWarnings(as.numeric(Y)),
      value = suppressWarnings(as.numeric(`Relative protein level`))
    )
  
  bad_x <- sum(!is.finite(df$x))
  bad_y <- sum(!is.finite(df$y))
  bad_value <- sum(!is.finite(df$value))
  
  if (bad_x > 0 || bad_y > 0 || bad_value > 0) {
    stop(
      paste0(
        "Non-finite values detected after reading CSV:\n",
        "bad x rows: ", bad_x, "\n",
        "bad y rows: ", bad_y, "\n",
        "bad value rows: ", bad_value, "\n",
        "请检查 CSV 是否有空单元格、文字、额外空格、合并单元格或非数字内容。"
      )
    )
  }
  
  return(df)
}

# =========================================================
# 2. Matrix bubble plot
# 2. 多行矩阵气泡图
# =========================================================
plot_xy_bubble <- function(
    df,
    
    # -------------------------
    # labels / 标签
    # -------------------------
    x_labels = NULL,   # 例如 c("NC_1","NC_2","KD_1","KD_2")
    y_labels = NULL,   # 例如 c("Protein1","Protein2","Protein3","Protein4","Protein5")
    
    # -------------------------
    # group annotation / 分组信息
    # -------------------------
    # 这三项最重要：
    # 1. group_names: 每个组显示什么名字
    # 2. group_xmin: 每组从哪一列开始
    # 3. group_xmax: 每组到哪一列结束
    #
    # 例如：
    # group_names = c("MM", "RE")
    # group_xmin  = c(1, 3)
    # group_xmax  = c(2, 4)
    #
    # 表示：
    # 第1组 "MM" 覆盖第1到第2列
    # 第2组 "RE" 覆盖第3到第4列
    #
    # 如果你有3组：
    # group_names = c("Ctrl", "MM", "RE")
    # group_xmin  = c(1, 3, 5)
    # group_xmax  = c(2, 4, 6)
    group_names = NULL,
    group_xmin = NULL,
    group_xmax = NULL,
    
    # -------------------------
    # plot settings / 作图参数
    # -------------------------
    plot_title = NULL,
    size_range = c(8, 26),
    legend_title = "Relative protein level",
    
    # 是否显示 x/y 标签
    show_x_labels = FALSE,
    show_y_labels = FALSE,
    show_group_line = TRUE,
    
    # -------------------------
    # style / 风格
    # -------------------------
    base_size = 18,
    plot_title_size = 22,
    axis_text_size = 16,
    legend_text_size = 16,
    legend_title_size = 18,
    
    # group title size / 分组标题字体大小
    group_text_size = 8,
    
    # 横线粗细
    group_line_width = 2.8,
    
    # 分组标题和横线的纵向位置偏移
    # 数值越大，标题越往上
    group_line_y_offset = 0.55,
    group_text_y_offset = 1.05,
    
    # 分组标题的横向整体偏移
    # 例如设为 0.2 可以整体向右挪
    group_text_x_offset = 0,
    
    # -------------------------
    # legend bar / 图例色条
    # -------------------------
    legend_bar_height_cm = 5.5,
    legend_bar_width_cm = 0.8,
    
    # -------------------------
    # margins / 边距
    # -------------------------
    top_margin = 25,
    right_margin = 20,
    bottom_margin = 20,
    left_margin = 20
) {
  
  if (any(!is.finite(df$x)) || any(!is.finite(df$y)) || any(!is.finite(df$value))) {
    stop("x / y / value contains non-finite values. 请先检查输入数据。")
  }
  
  x_vals <- sort(unique(df$x))
  y_vals <- sort(unique(df$y))
  
  # 让 Y=1 在最上面
  max_y <- max(y_vals)
  df2 <- df %>%
    mutate(y_plot = max_y - y + 1)
  
  if (is.null(x_labels)) x_labels <- as.character(x_vals)
  if (is.null(y_labels)) y_labels <- as.character(y_vals)
  
  y_labels_plot <- rev(y_labels)
  y_plot_vals <- sort(unique(df2$y_plot))
  
  p <- ggplot(df2, aes(x = x, y = y_plot)) +
    geom_point(
      aes(size = value, fill = value),
      shape = 21,
      stroke = 0
    ) +
    scale_size_continuous(range = size_range) +
    
    # 蓝色马卡龙大跨度渐变
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#DCEBFA", "#AFCFEF", "#6FA8DC", "#2F6FA3"),
      name = legend_title
    ) +
    
    guides(
      size = "none",
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0,
        barheight = unit(legend_bar_height_cm, "cm"),
        barwidth = unit(legend_bar_width_cm, "cm"),
        frame.colour = NA,
        ticks.colour = "white"
      )
    ) +
    labs(title = plot_title) +
    theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = plot_title_size,
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "right",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      plot.margin = margin(top_margin, right_margin, bottom_margin, left_margin)
    )
  
  # -------------------------
  # group lines and titles
  # 分组横线和标题
  # -------------------------
  if (show_group_line && !is.null(group_names) && !is.null(group_xmin) && !is.null(group_xmax)) {
    
    # 检查三者长度是否一致
    if (!(length(group_names) == length(group_xmin) && length(group_names) == length(group_xmax))) {
      stop("group_names, group_xmin, group_xmax 长度必须一致。")
    }
    
    group_df <- data.frame(
      group = group_names,
      xmin = group_xmin,
      xmax = group_xmax
    ) %>%
      mutate(
        xmid = (xmin + xmax) / 2 + group_text_x_offset,
        y_line = max(df2$y_plot) + group_line_y_offset,
        y_text = max(df2$y_plot) + group_text_y_offset
      )
    
    p <- p +
      geom_segment(
        data = group_df,
        aes(
          x = xmin - 0.15,
          xend = xmax + 0.15,
          y = y_line,
          yend = y_line
        ),
        inherit.aes = FALSE,
        linewidth = group_line_width,
        color = "black"
      ) +
      geom_text(
        data = group_df,
        aes(
          x = xmid,
          y = y_text,
          label = group
        ),
        inherit.aes = FALSE,
        size = group_text_size
      )
  }
  
  # -------------------------
  # axis-like labels
  # 类坐标轴标签
  # -------------------------
  p <- p +
    scale_x_continuous(
      breaks = if (show_x_labels) x_vals else NULL,
      labels = if (show_x_labels) x_labels else NULL,
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    scale_y_continuous(
      breaks = if (show_y_labels) y_plot_vals else NULL,
      labels = if (show_y_labels) y_labels_plot else NULL,
      expand = expansion(mult = c(0.03, 0.22))
    ) +
    coord_cartesian(clip = "off")
  
  if (show_x_labels || show_y_labels) {
    p <- p +
      theme(
        axis.text.x = if (show_x_labels) {
          element_text(size = axis_text_size, color = "black")
        } else {
          element_blank()
        },
        axis.text.y = if (show_y_labels) {
          element_text(size = axis_text_size, color = "black")
        } else {
          element_blank()
        }
      )
  }
  
  return(p)
}

# =========================================================
# 3. Single-row bubble plot
# 3. 单行气泡图（适合只有 y=1 的数据）
# =========================================================
plot_xy_bubble_row <- function(
    df,
    x_labels = NULL,
    plot_title = NULL,
    legend_title = NULL,
    big_bottom_title = "Relative Intensity(AU)",
    size_range = c(8, 28),
    show_top_axis = TRUE,
    show_left_colorbar = TRUE,
    base_size = 16,
    big_title_size = 34,
    legend_text_size = 14,
    legend_title_size = 16,
    top_margin = 20,
    right_margin = 20,
    bottom_margin = 20,
    left_margin = 20
) {
  
  # 只保留第一种 y
  y_first <- sort(unique(df$y))[1]
  df2 <- df %>%
    filter(y == y_first)
  
  x_vals <- sort(unique(df2$x))
  if (is.null(x_labels)) x_labels <- as.character(x_vals)
  
  df2 <- df2 %>%
    mutate(y_plot = 1)
  
  p <- ggplot(df2, aes(x = x, y = y_plot)) +
    geom_point(
      aes(size = value, fill = value),
      shape = 21,
      stroke = 0
    ) +
    scale_size_continuous(range = size_range) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#DCEBFA", "#AFCFEF", "#6FA8DC", "#2F6FA3"),
      name = legend_title %||% ""
    ) +
    guides(
      size = "none",
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0,
        barheight = unit(3.5, "cm"),
        barwidth = unit(0.8, "cm"),
        frame.colour = NA,
        ticks.colour = "white"
      )
    ) +
    labs(title = plot_title) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = if (show_left_colorbar) "left" else "none",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      plot.margin = margin(top_margin, right_margin, bottom_margin, left_margin)
    ) +
    scale_x_continuous(
      breaks = x_vals,
      labels = x_labels,
      expand = expansion(mult = c(0.03, 0.03)),
      position = if (show_top_axis) "top" else "bottom"
    ) +
    coord_cartesian(clip = "off") +
    annotate(
      "text",
      x = mean(range(df2$x)),
      y = 0.35,
      label = big_bottom_title,
      size = big_title_size / 4
    )
  
  if (show_top_axis) {
    p <- p +
      theme(
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_line(color = "black", linewidth = 0.8),
        axis.line.x.top = element_line(color = "black", linewidth = 0.8)
      )
  }
  
  return(p)
}

# =========================================================
# 4. Auto mode
# 4. 自动判断模式
# =========================================================
# 如果只有一个 y 值，就自动画单行图
# 如果有多个 y 值，就画矩阵图
plot_xy_bubble_auto <- function(
    df,
    x_labels = NULL,
    y_labels = NULL,
    group_names = NULL,
    group_xmin = NULL,
    group_xmax = NULL,
    legend_title = "Relative protein level",
    row_bottom_title = "Relative Intensity(AU)"
) {
  n_y <- length(unique(df$y))
  
  if (n_y == 1) {
    p <- plot_xy_bubble_row(
      df,
      x_labels = x_labels,
      legend_title = legend_title,
      big_bottom_title = row_bottom_title
    )
  } else {
    p <- plot_xy_bubble(
      df,
      x_labels = x_labels,
      y_labels = y_labels,
      group_names = group_names,
      group_xmin = group_xmin,
      group_xmax = group_xmax,
      legend_title = legend_title,
      show_group_line = TRUE
    )
  }
  
  return(p)
}

# =========================================================
# 5. Save helper
# 5. 保存函数
# =========================================================
save_bubble_plot <- function(
    plot_obj,
    output_prefix = "bubble_plot",
    width = 8,
    height = 6,
    dpi = 600
) {
  ggsave(
    filename = paste0(output_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(output_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
}

# =========================================================
# 6. Example usage
# 6. 使用示例
# =========================================================
bubble_df <- read_xy_bubble_data("western_blot_xy.csv")

print(bubble_df)
summary(bubble_df)

# -------------------------
# Example A: matrix mode
# 示例A：矩阵模式
# -------------------------
p1 <- plot_xy_bubble(
  bubble_df,
  
  # x轴每列名字
  x_labels = c("NC_1", "NC_2", "KD_1", "KD_2"),
  
  # y轴每行名字
  y_labels = c("Protein1", "Protein2", "Protein3", "Protein4", "Protein5"),
  
  # -------------------------
  # group 设置
  # -------------------------
  # 现在是 2 组：
  # 第1组叫 MM，覆盖第1-2列
  # 第2组叫 RE，覆盖第3-4列
  group_names = c("MM", "RE"),
  group_xmin = c(1, 3),
  group_xmax = c(2, 4),
  
  # 如果你想改成 3 组，例如：
  # group_names = c("Ctrl", "MM", "RE")
  # group_xmin = c(1, 3, 5)
  # group_xmax = c(2, 4, 6)
  
  plot_title = NULL,
  legend_title = "Relative protein level",
  show_x_labels = FALSE,
  show_y_labels = FALSE,
  show_group_line = TRUE,
  
  # 分组标题位置和样式
  group_text_size = 8,
  group_line_width = 2.8,
  group_line_y_offset = 0.55,  # 横线往上/往下
  group_text_y_offset = 1.05,  # 标题往上/往下
  group_text_x_offset = 0      # 标题整体左右移动
)

print(p1)

save_bubble_plot(
  p1,
  output_prefix = "protein_bubble_xy_blue",
  width = 9,
  height = 7,
  dpi = 600
)

