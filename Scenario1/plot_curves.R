library(ggplot2)
library(dplyr)
library(tidyr)
load("plot_curves_only.RData")

h <- 0.01
designpoints <- seq(0, 1, by = h)

n_func      <- 10
func_labels <- sprintf("italic(l)==%d", seq_len(n_func))

mat_list <- list(
  "Truth"   = ISE_true,
  "Lasso"   = est_coef_list[[1]],
  "SCAD"    = est_coef_list[[2]],
  "MCP"     = est_coef_list[[3]],
  "FLCRM"   = est_coef_list[[4]],
  "Lasso-w" = est_coef_list[[5]],
  "SCAD-w"  = est_coef_list[[6]],
  "MCP-w"   = est_coef_list[[7]]
)

to_df <- function(m, method) {
  m <- as.matrix(m)
  m <- m[, seq_len(n_func), drop = FALSE]
  df <- data.frame(m, check.names = FALSE)
  names(df) <- func_labels
  df$x <- designpoints
  df %>%
    pivot_longer(cols = all_of(func_labels),
                 names_to = "func", values_to = "y") %>%
    mutate(Method = method)
}

df_all <- bind_rows(lapply(names(mat_list),
                           function(m) to_df(mat_list[[m]], m)))

df_all$func   <- factor(df_all$func,   levels = func_labels)
df_all$Method <- factor(df_all$Method,
                        levels = c( "Lasso", "Lasso-w",
                                    "SCAD", "SCAD-w",
                                    "MCP", "MCP-w",
                                    "Truth","FLCRM"))

my_colors <- c(
  "Truth"   = "black",
  "Lasso"   = "#E41A1C",
  "Lasso-w" = "#FC9272",
  "SCAD"    = "#377EB8",
  "SCAD-w"  = "#9ECAE1",
  "MCP"     = "#4DAF4A",
  "MCP-w"   = "#A1D99B",
  "FLCRM"   = "#984EA3"
)

my_linetypes <- c(
  "Truth"   = "solid",
  "Lasso"   = "dashed",
  "Lasso-w" = "dashed",
  "SCAD"    = "dotdash",
  "SCAD-w"  = "dotdash",
  "MCP"     = "dotted",
  "MCP-w"   = "dotted",
  "FLCRM"   = "longdash"
)

my_sizes <- c(
  "Truth" = 0.9, "Lasso" = 0.6, "Lasso-w" = 0.6,
  "SCAD"  = 0.6, "SCAD-w" = 0.6,
  "MCP"   = 0.6, "MCP-w"  = 0.6,
  "FLCRM" = 0.6
)

p <- ggplot(df_all,
            aes(x = x, y = y,
                color = Method, linetype = Method, linewidth = Method)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ func, ncol = 2, scales = "free_y",
             labeller = labeller(func = label_parsed)) +
  scale_color_manual(values = my_colors,       name = "Method") +
  scale_linetype_manual(values = my_linetypes, name = "Method") +
  scale_linewidth_manual(values = my_sizes,    name = "Method") +
  labs(x = NULL, y = " ") +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

ann <- data.frame(
  func = factor(func_labels[c(9, 10)], levels = func_labels),
  x    = 0.5,
  y    = -Inf,
  lab  = "italic(s)"
)

p <- p +
  geom_text(data = ann,
            aes(x = x, y = y, label = lab),
            inherit.aes = FALSE,
            parse    = TRUE,
            vjust    = 3,
            size     = 4.2,
            color    = "black") +
  coord_cartesian(clip = "off")

print(p)





