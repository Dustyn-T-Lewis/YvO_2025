################################################################################
#   Figure 0 — Legend Keys
#   Standalone legends for F0 phenotype panels (training volume, body comp)
################################################################################

source("04_Figures_v2/keys/extract_keys.R")

message("Building F0 keys...")

# --- Group fill legend (Young_Pre, Young_Post, Old_Pre, Old_Post) ---
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

# --- Assemble ---
keys_F0 <- wrap_elements(key_group)

save_key(keys_F0, "F0_keys", width = 280, height = 60)

message("  F0 keys done")
