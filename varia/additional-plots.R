library(tidyverse)
library(ggforce)

# Prentice Definition -------------------------------------------------------

# Function that will demarcate the region in which trial-level treatment effects
# are allowed to live to satisfy Prentice's definition.
upper_f = function(x) {
  ifelse(x >= 0, 100 * x ^ 2, -sqrt(-0.1 * x))
}

lower_f = function(x) {
  ifelse(x < 0, -100 * x ^ 2, +sqrt(0.1 * x))
}

x = seq(from = -3,
        to = 3,
        length.out = 1e3)
plotting_data =
  tibble(x, y_upper = upper_f(x), y_lower = lower_f(x))

ggplot(plotting_data, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = y_upper), size = 1) +
  geom_line(aes(y = y_lower), size = 1) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(x = 0, y = 0), size = 3, color = "blue") +
  labs(
    title = "Prentice's Definition",
    x = expr(alpha), y = expr(beta)
  )

ggsave(
  filename = "prentice.pdf",
  device = "pdf",
  path = "varia",
  width = single_width,
  height = single_height, 
  units = unit
)

# Trial-level Surrogacy ------------------------------------------------------

# Define a sequence of ellipses
ellipse_data <- tibble(
  x0 = 0,      # ellipse centers (x)
  y0 = 0,                                # all on horizontal line y = 0
  size = seq(1, 2.5, length.out = 3),   # size factor
  a = size,                              # major axis
  b = size * 0.15,                        # minor axis (same ratio)
  angle = pi / 4                            # no rotation (horizontal major axis)
)

# Plot
ggplot() +
  geom_ellipse(data = ellipse_data,
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle),
               fill = "lightblue", alpha = 0.5, color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 2)) +
  labs(
    title = "Trial-Level Surrogacy",
    x = expr(alpha), y = expr(beta)
  )

ggsave(
  filename = "trial-level-surrogacy.pdf",
  device = "pdf",
  path = "varia",
  width = single_width,
  height = single_height, 
  units = unit
)
