#封装函数
plot_single_protein_bubble_clean <- function(
    df,
    bottom_title = "Relative Intensity(AU)",
    legend_title = NULL,
    
    # 蓝色马卡龙大跨度渐变
    bubble_colours = c("#F7FBFF", "#DCEBFA", "#AFCFEF", "#6FA8DC", "#2F6FA3"),
    
    # 点大小
    size_range = c(8, 30),
    
    # 自定义 x 坐标，可压缩点间距
    custom_x = NULL,
    
    # 布局
    bubble_y = 0.72,
    axis_y = 0.84,
    title_y = 0.34,
    
    # 顶部轴
    show_top_axis = TRUE,
    major_breaks = NULL,
    minor_breaks = NULL,
    major_tick_len = 0.035,
    minor_tick_len = 0.020,
    axis_line_width = 0.9,
    tick_line_width = 0.8,
    
    # 字体
    bottom_title_size = 34,
    legend_text_size = 16,
    legend_title_size = 16,
    
    # 图例
    show_left_colorbar = TRUE,
    legend_bar_height_cm = 2.8,
    legend_bar_width_cm = 0.55,
    
    # 图例与画图区的距离
    legend_box_margin_pt = 6,
    
    # 整体范围
    x_expand = 0.10,
    y_limits = c(0.18, 0.92),
    
    # 页边距
    plot_margin_top = 10,
    plot_margin_right = 20,
    plot_margin_bottom = 10,
    plot_margin_left = 30
) {
  
  first_y <- sort(unique(df$y))[1]
  df2 <- df %>%
    filter(y == first_y) %>%
    arrange(x)
  
  # 可选：自定义横坐标，用来压缩点间距
  if (!is.null(custom_x)) {
    if (length(custom_x) != nrow(df2)) {
      stop("custom_x 的长度必须等于单行点的数量。")
    }
    df2$x_plot <- custom_x
  } else {
    df2$x_plot <- df2$x
  }
  
  df2$y_plot <- bubble_y
  
  x_vals <- df2$x_plot
  
  if (is.null(major_breaks)) {
    major_breaks <- x_vals
  }
  if (is.null(minor_breaks)) {
    minor_breaks <- numeric(0)
  }
  
  axis_df <- data.frame(
    x = min(x_vals) - 0.10,
    xend = max(x_vals) + 0.10,
    y = axis_y,
    yend = axis_y
  )
  
  if (length(major_breaks) > 0) {
    major_ticks_df <- data.frame(
      x = major_breaks,
      xend = major_breaks,
      y = rep(axis_y, length(major_breaks)),
      yend = rep(axis_y + major_tick_len, length(major_breaks))
    )
  } else {
    major_ticks_df <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
  }
  
  if (length(minor_breaks) > 0) {
    minor_ticks_df <- data.frame(
      x = minor_breaks,
      xend = minor_breaks,
      y = rep(axis_y, length(minor_breaks)),
      yend = rep(axis_y + minor_tick_len, length(minor_breaks))
    )
  } else {
    minor_ticks_df <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
  }
  
  p <- ggplot(df2, aes(x = x_plot, y = y_plot)) +
    geom_point(
      aes(size = value, fill = value),
      shape = 21,
      stroke = 0
    ) +
    scale_size_continuous(range = size_range) +
    scale_fill_gradientn(
      colours = bubble_colours,
      name = legend_title %||% ""
    ) +
    guides(
      size = "none",
      fill = if (show_left_colorbar) {
        guide_colorbar(
          title.position = "top",
          title.hjust = 0,
          barheight = unit(legend_bar_height_cm, "cm"),
          barwidth = unit(legend_bar_width_cm, "cm"),
          frame.colour = NA,
          ticks.colour = "white"
        )
      } else {
        "none"
      }
    ) +
    
    # 顶部横线
    {if (show_top_axis) geom_segment(
      data = axis_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = axis_line_width,
      color = "black"
    )} +
    
    # 主刻度
    {if (show_top_axis && nrow(major_ticks_df) > 0) geom_segment(
      data = major_ticks_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = tick_line_width,
      color = "black"
    )} +
    
    # 次刻度
    {if (show_top_axis && nrow(minor_ticks_df) > 0) geom_segment(
      data = minor_ticks_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = tick_line_width,
      color = "black"
    )} +
    
    # 下方大标题
    annotate(
      "text",
      x = mean(range(x_vals)),
      y = title_y,
      label = bottom_title,
      size = bottom_title_size / 4
    ) +
    
    scale_x_continuous(
      limits = c(min(x_vals) - x_expand, max(x_vals) + x_expand),
      breaks = NULL,
      labels = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = NULL,
      labels = NULL,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 16) +
    theme(
      legend.position = if (show_left_colorbar) "left" else "none",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.box.margin = margin(0, legend_box_margin_pt, 0, 0),
      plot.margin = margin(
        plot_margin_top,
        plot_margin_right,
        plot_margin_bottom,
        plot_margin_left
      )
    )
  
  return(p)
}

##调用函数

single_df <- read_single_protein_data("singel protein.csv")

p_single <- plot_single_protein_bubble_clean(
  single_df,
  bottom_title = "Relative Intensity(AU)",
  legend_title = NULL,
  size_range = c(15, 30),
  
  # 压缩点间距
  custom_x = c(1, 1.9, 2.8, 3.7),
  
  bubble_y = 0.65,
  axis_y = 0.72,
  title_y = 0.47,
  
  # 每个点一个长刻度
  major_breaks = c(1, 1.9, 2.8, 3.7),
  minor_breaks = NULL,
  
  show_left_colorbar = TRUE,
  legend_bar_height_cm = 3.0,
  legend_bar_width_cm = 0.8,
  legend_box_margin_pt = 10,
  
  bottom_title_size = 34,
  x_expand = 0.10,
  y_limits = c(0.18, 0.92),
  
  # 给左边图例更多空间
  plot_margin_left = 60
)

print(p_single)


###保存函数##
save_single_protein_plot <- function(
    plot_obj,
    output_prefix = "single_protein_bubble",
    width = 8,
    height = 4,
    dpi = 600
) {
  png_file <- paste0(output_prefix, ".png")
  pdf_file <- paste0(output_prefix, ".pdf")
  
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  message("Saved PNG: ", normalizePath(png_file, mustWork = FALSE))
  
  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  message("Saved PDF: ", normalizePath(pdf_file, mustWork = FALSE))
}
##保存##
save_single_protein_plot(
  p_single,
  output_prefix = "single_protein_bubble_final",
  width = 8,
  height = 6,
  dpi = 600
)
