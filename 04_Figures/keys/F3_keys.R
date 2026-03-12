# Figure 3 — Legend Keys
source("04_Figures/keys/extract_keys.R")

# Direction colors (volcano rings)
dir_df <- tibble(
  direction = factor(c("Up", "Down", "NS"), levels = c("Up", "Down", "NS")),
  y = 1:3
)
p_dir <- ggplot(dir_df, aes(x = 1, y = y, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = DIR_COLORS, name = "Direction") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_dir <- cowplot::get_legend(p_dir)

# Reversal quadrant colors
rev_df <- tibble(
  quadrant = factor(c("Rejuvenation", "Reversal",
                       "Non-restored (+/+)", "Non-restored (-/-)"),
                    levels = c("Rejuvenation", "Reversal",
                               "Non-restored (+/+)", "Non-restored (-/-)")),
  y = 1:4
)
rev_colors <- c("Rejuvenation" = "#4393C3", "Reversal" = "#D6604D",
                 "Non-restored (+/+)" = "#E0E0E0", "Non-restored (-/-)" = "#B0B0B0")
p_rev <- ggplot(rev_df, aes(x = 1, y = y, fill = quadrant)) +
  geom_col() +
  scale_fill_manual(values = rev_colors, name = "Reversal Category") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 9, face = "bold"))
key_rev <- cowplot::get_legend(p_rev)

# RRHO2 -log10(p) colorbar
rrho_df <- expand.grid(x = 1:10, y = 1:10)
rrho_df$value <- runif(100, 0, 50)
p_rrho <- ggplot(rrho_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "#FFFFCC", "#FD8D3C", "#E31A1C", "#800026"),
    name = expression(-log[10](p)),
    limits = c(0, 50)
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.width = unit(20, "mm"),
        legend.key.height = unit(3, "mm"))
key_rrho <- cowplot::get_legend(p_rrho)

keys_F3 <- wrap_elements(key_dir) + wrap_elements(key_rev) +
           wrap_elements(key_rrho) +
  plot_layout(nrow = 1)

save_key(keys_F3, "F3_keys", width = 350, height = 70)
