################################################################################
#   Figure 6 — Legend Keys
#   Standalone legends for F6 translational utility panels
################################################################################

source("04_Figures_v2/keys/extract_keys.R")

message("Building F6 keys...")

# --- Age colors (PCA, scatter panels B/C) ---
age_df <- tibble(
  age = factor(c("Young", "Old"), levels = c("Young", "Old")),
  y = 1:2
)
p_age <- ggplot(age_df, aes(x = 1, y = y, color = age, fill = age)) +
  geom_point(size = 3) +
  scale_color_manual(values = AGE_COLORS, name = "Age") +
  scale_fill_manual(values = AGE_COLORS, name = "Age") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_age <- cowplot::get_legend(p_age)

# --- Pearson r colorbar (heatmap, Panel E) ---
heat_df <- expand.grid(x = 1:10, y = 1)
heat_df$r <- seq(-1, 1, length.out = 10)
p_heat <- ggplot(heat_df, aes(x = x, y = y, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.width = unit(25, "mm"),
        legend.key.height = unit(3, "mm"))
key_heat <- cowplot::get_legend(p_heat)

# --- Timepoint legend (ROC comparison) ---
tp_roc_df <- tibble(
  timepoint = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
  x = 1:2, y = 1:2
)
p_tp <- ggplot(tp_roc_df, aes(x = x, y = y, color = timepoint, linetype = timepoint)) +
  geom_line() +
  scale_color_manual(values = c(Pre = "#D6604D", Post = "#5DA5DA"),
                      name = "Timepoint") +
  scale_linetype_manual(values = c(Pre = "solid", Post = "dashed"),
                         name = "Timepoint") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_tp <- cowplot::get_legend(p_tp)

# --- Significance notation ---
sig_df <- tibble(
  label = c("* p_adj < .05", "** p_adj < .01", "*** p_adj < .001"),
  x = 1:3, y = 1
)
p_sig <- ggplot(sig_df, aes(x = x, y = y, label = label)) +
  geom_text(size = 3, fontface = "bold") +
  theme_void() +
  labs(title = "BH-corrected significance") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# --- Assemble ---
keys_F6 <- wrap_elements(key_age) + wrap_elements(key_heat) +
           wrap_elements(key_tp) + wrap_elements(p_sig) +
  plot_layout(nrow = 1, widths = c(1, 2, 1, 1.5))

save_key(keys_F6, "F6_keys", width = 350, height = 70)

message("  F6 keys done")
