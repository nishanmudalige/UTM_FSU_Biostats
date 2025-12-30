library(ggplot2)
library(dplyr)
library(purrr)
library(grid)

set.seed(1)

# ---- tuned to match the look ----
n_paths <- 10
y0      <- 0.16
x_min   <- 0
x_max   <- 1.08
y_min   <- 0
y_max   <- 1.12

# Approx palette similar to the image
pal <- c(
  "#F43F5E",
  "#FB7185",
  "#F97316",
  "#FBBF24",
  "#EAB308",
  "#22C55E",
  "#14B8A6",
  "#06B6D4",
  "#2563EB",
  "#7C3AED"
)

x <- seq(0, 1, length.out = 420)

# helper: smooth random noise (no extra packages)
smooth_noise <- function(n, sd = 1, k = 41) {
  z <- rnorm(n, sd = sd)
  t <- seq(-2, 2, length.out = k)
  ker <- dnorm(t)
  ker <- ker / sum(ker)
  as.numeric(stats::filter(z, ker, sides = 2, circular = FALSE))
}

make_path <- function(id) {
  
  y_end <- sample(
    c(runif(1, 0.05, 0.25),
      runif(1, 0.25, 0.55),
      runif(1, 0.55, 1.05)),
    size = 1,
    prob = c(0.18, 0.42, 0.40)
  )
  
  baseline <- y0 + (y_end - y0) * x
  
  env <- (1 - exp(-6 * x)) * (0.85 + 0.15 * (1 - x))
  
  m     <- sample(3:5, 1)
  freqs <- sample(1:4, m, replace = TRUE)
  amps  <- runif(m, 0.03, 0.10) * sample(c(0.8, 1, 1.2), m, TRUE)
  phase <- runif(m, 0, 2*pi)
  
  wiggle <- rowSums(sapply(seq_len(m), function(j) {
    amps[j] * sin(2*pi*freqs[j]*x + phase[j])
  }))
  
  nse <- smooth_noise(length(x), sd = runif(1, 0.03, 0.06), k = 51)
  nse[is.na(nse)] <- 0
  
  y <- baseline + env * (wiggle + 0.55 * nse)
  
  y[1] <- y0
  y <- pmin(pmax(y, 0.02), 1.08)
  
  tibble(id = factor(id), x = x, y = y)
}

df <- map_dfr(1:n_paths, make_path)

p <- ggplot(df, aes(x, y, group = id, colour = id)) +
  geom_line(linewidth = 2.6, lineend = "round") +
  geom_point(
    data = tibble(x = 0, y = y0),
    aes(x, y),
    inherit.aes = FALSE,
    colour = "#1D4ED8",
    size = 3.8
  ) +
  annotate(
    "segment",
    x = 0, y = 0, xend = x_max, yend = 0,
    linewidth = 3.2, colour = "#C7C7C7", lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.30, "in"), angle = 25)
  ) +
  annotate(
    "segment",
    x = 0, y = 0, xend = 0, yend = y_max,
    linewidth = 3.2, colour = "#C7C7C7", lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.30, "in"), angle = 25)
  ) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = seq(0, 1, by = 0.10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(0, 1, by = 0.10),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    # ---- TRANSPARENT BACKGROUND ----
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    
    panel.grid.major = element_line(colour = "#EEF0F3", linewidth = 0.9),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(40, 55, 45, 55),
    aspect.ratio = 0.62
  )

p

# ---- SAVE with transparency (recommended) ----
# install.packages("ragg")  # if needed
# ragg::agg_png("plot_transparent.png", width = 1600, height = 1000, res = 200, background = "transparent")
res = 1200
ragg::agg_png("plot_transparent.png", 
              width = 1600*(res/200), 
              height = 1000*(res/200), res = res, background = "transparent")
print(p)
dev.off()

# Alternative (also works on most systems):
# ggsave("plot_transparent.png", p, bg = "transparent", width = 8, height = 5, dpi = 300)
