# Figure 5 — Legend Keys
source("04_Figures/keys/extract_keys.R")

# Pearson r colorbar (heatmap, Panel B)
heat_df <- expand.grid(x = 1:10, y = 1)
heat_df$r <- seq(-0.8, 0.8, length.out = 10)
p_heat <- ggplot(heat_df, aes(x = x, y = y, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0033CC", mid = "white", high = "#CC0000",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       name = "Pearson r") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.width = unit(25, "mm"),
        legend.key.height = unit(3, "mm"))
key_heat <- cowplot::get_legend(p_heat)

# Significance notation
sig_df <- tibble(
  label = c("* p < .05", "** p < .01", "*** p < .001"),
  x = 1:3, y = 1
)
p_sig <- ggplot(sig_df, aes(x = x, y = y, label = label)) +
  geom_text(size = 3, fontface = "bold") +
  theme_void() +
  labs(title = "Significance") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# Module preservation thresholds
pres_df <- tibble(
  label = c("Zsummary > 10: Strong", "Zsummary 2-10: Moderate", "Zsummary < 2: Not preserved"),
  x = 1, y = 3:1
)
p_pres <- ggplot(pres_df, aes(x = x, y = y, label = label)) +
  geom_text(size = 2.5, hjust = 0) +
  theme_void() +
  labs(title = "Preservation") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# Age colors (triptych scatter)
age_df <- tibble(
  age = factor(c("Young", "Old"), levels = c("Young", "Old")),
  y = 1:2
)
p_age <- ggplot(age_df, aes(x = 1, y = y, color = age)) +
  geom_point(size = 3) +
  scale_color_manual(values = AGE_COLORS, name = "Age") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_age <- cowplot::get_legend(p_age)

keys_F5 <- wrap_elements(key_heat) + wrap_elements(key_age) +
           wrap_elements(p_sig) + wrap_elements(p_pres) +
  plot_layout(nrow = 1, widths = c(2, 1, 1.5, 2))

save_key(keys_F5, "F5_keys", width = 350, height = 70)
