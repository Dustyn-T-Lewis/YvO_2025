# Figure 1 — Legend Keys
source("04_Figures/keys/extract_keys.R")

# Group fill (CV violins, Panel A)
group_df <- tibble(
  group = factor(names(GROUP_FILL), levels = names(GROUP_FILL)),
  y = seq_along(GROUP_FILL)
)
p_group <- ggplot(group_df, aes(x = 1, y = y, fill = group)) +
  geom_col() +
  scale_fill_manual(values = GROUP_FILL, name = "Group") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_group <- cowplot::get_legend(p_group)

# Contrast colors (logFC density, Panel B)
ctr_df <- tibble(
  contrast = factor(c("Aging", "Training_Young", "Training_Old"),
                    levels = c("Aging", "Training_Young", "Training_Old")),
  y = 1:3
)
p_contrast <- ggplot(ctr_df, aes(x = 1, y = y, fill = contrast)) +
  geom_col() +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")],
                    name = "Contrast") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_contrast <- cowplot::get_legend(p_contrast)

# PCA color/shape (Panel C)
pca_df <- tibble(
  group = factor(names(PCA_COLORS), levels = names(PCA_COLORS)),
  x = 1:4, y = 1:4
)
p_pca <- ggplot(pca_df, aes(x = x, y = y, color = group, shape = group, fill = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = PCA_COLORS, name = "Group") +
  scale_fill_manual(values = PCA_COLORS, name = "Group") +
  scale_shape_manual(values = PCA_SHAPES, name = "Group") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_pca <- cowplot::get_legend(p_pca)

# Direction colors (fGSEA + UpSet, Panels E/F)
dir_df <- tibble(
  direction = factor(c("Up", "Down"), levels = c("Up", "Down")),
  y = 1:2
)
p_dir <- ggplot(dir_df, aes(x = 1, y = y, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = DIR_COLORS[c("Up", "Down")], name = "Direction") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_dir <- cowplot::get_legend(p_dir)

keys_F1 <- wrap_elements(key_group) + wrap_elements(key_contrast) +
           wrap_elements(key_pca) + wrap_elements(key_dir) +
  plot_layout(nrow = 1)

save_key(keys_F1, "F1_keys", width = 350, height = 60)
