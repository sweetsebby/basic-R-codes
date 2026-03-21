suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

dir.create("results/ci", recursive = TRUE, showWarnings = FALSE)

# 优先尝试读取你仓库里现有的 western_blot_xy.csv
# 如果没有，或者格式不合适，就自动生成一份 toy data
csv_file <- "western_blot_xy.csv"

make_toy_data <- function() {
  data.frame(
    x = seq(1, 50),
    y = cumsum(rnorm(50, mean = 0.05, sd = 0.2))
  )
}

if (file.exists(csv_file)) {
  dat <- tryCatch(
    read_csv(csv_file, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (!is.null(dat)) {
    numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
    if (length(numeric_cols) >= 2) {
      plot_dat <- data.frame(
        x = dat[[numeric_cols[1]]],
        y = dat[[numeric_cols[2]]]
      )
    } else {
      plot_dat <- make_toy_data()
    }
  } else {
    plot_dat <- make_toy_data()
  }
} else {
  plot_dat <- make_toy_data()
}

p <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_line(linewidth = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = "CI example figure",
    x = "X",
    y = "Y"
  )

ggsave(
  filename = "results/ci/ci_example_figure.pdf",
  plot = p,
  width = 4,
  height = 3
)

ggsave(
  filename = "results/ci/ci_example_figure.png",
  plot = p,
  width = 4,
  height = 3,
  dpi = 300
)

stopifnot(file.exists("results/ci/ci_example_figure.pdf"))
stopifnot(file.exists("results/ci/ci_example_figure.png"))

message("CI figure rendering finished successfully.")
