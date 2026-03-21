suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# --------- CI example output directory ---------
out_dir <- file.path("results", "ci")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------- where to read example CSV ---------
csv_file <- file.path("data", "example", "western_blot_xy.csv")
stopifnot(file.exists(csv_file))

# 优先读取示例 CSV；读失败就用 toy data
dat <- tryCatch(
  read_csv(csv_file, show_col_types = FALSE),
  error = function(e) NULL
)

make_toy_data <- function() {
  data.frame(
    x = seq(1, 50),
    y = cumsum(rnorm(50, mean = 0.05, sd = 0.2))
  )
}

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

p <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_line(linewidth = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = "CI example figure",
    x = "X",
    y = "Y"
  )

ggsave(
  filename = file.path(out_dir, "ci_example_figure.pdf"),
  plot = p,
  width = 4,
  height = 3
)

ggsave(
  filename = file.path(out_dir, "ci_example_figure.png"),
  plot = p,
  width = 4,
  height = 3,
  dpi = 300
)

stopifnot(file.exists(file.path(out_dir, "ci_example_figure.pdf")))
stopifnot(file.exists(file.path(out_dir, "ci_example_figure.png")))

message("CI figure rendering finished successfully.")
