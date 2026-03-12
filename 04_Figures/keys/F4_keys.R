# Figure 4 — Legend Keys
source("04_Figures/keys/extract_keys.R")

# FCM cluster colors
fcm_pal <- c("C1" = "#F8766D", "C2" = "#7CAE00", "C3" = "#00BFC4", "C4" = "#C77CFF")
fcm_df <- tibble(
  cluster = factor(names(fcm_pal), levels = names(fcm_pal)),
  y = seq_along(fcm_pal)
)
p_fcm <- ggplot(fcm_df, aes(x = 1, y = y, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = fcm_pal, name = "FCM Cluster") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_fcm <- cowplot::get_legend(p_fcm)

# Age colors (PCA scatter)
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

# Timepoint shapes
tp_df <- tibble(
  time = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
  x = 1:2, y = 1:2
)
p_tp <- ggplot(tp_df, aes(x = x, y = y, shape = time)) +
  geom_point(size = 3) +
  scale_shape_manual(values = SHAPE_TP, name = "Timepoint") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
key_tp <- cowplot::get_legend(p_tp)

keys_F4 <- wrap_elements(key_fcm) + wrap_elements(key_age) +
           wrap_elements(key_tp) +
  plot_layout(nrow = 1)

save_key(keys_F4, "F4_keys", width = 300, height = 60)
