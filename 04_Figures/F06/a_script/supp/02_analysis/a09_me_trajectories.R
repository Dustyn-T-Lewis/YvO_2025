# F06 SUPP — Module eigengene trajectories (Pre→Post, faceted by Age)
#
# Complements Panel A's ΔYoung / ΔOld columns by showing the underlying
# per-subject Pre→Post ME trajectories for each non-grey module.
#
# Output: 04_Figures/F06/b_reports/supp/02_analysis/SUPP_me_trajectories.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)

DAT <- "04_Figures/F06/c_data"
RPT <- "04_Figures/F06/b_reports/supp/02_analysis"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

MEs  <- readRDS(file.path(DAT, "MEs.rds"))
meta <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE) %>%
  mutate(time = factor(time, levels = c("Pre", "Post")),
         age  = factor(age,  levels = c("Young", "Old")))
mod_bio_labels <- read_csv(file.path(DAT, "mod_bio_labels.csv"),
                           show_col_types = FALSE)

non_grey <- setdiff(colnames(MEs), "MEgrey")
mod_color_raw <- setNames(sub("^ME", "", non_grey), non_grey)
bio_map       <- setNames(mod_bio_labels$bio_label,
                          paste0("ME", mod_bio_labels$module_color))
display_label <- ifelse(!is.na(bio_map[non_grey]),
                        sprintf("%s\n(%s)",
                                bio_map[non_grey],
                                str_to_title(mod_color_raw)),
                        str_to_title(mod_color_raw))
names(display_label) <- non_grey

# Long-format trajectories: one row per (subject, timepoint, module)
traj <- MEs %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(all_of(non_grey), names_to = "module", values_to = "ME") %>%
  left_join(meta %>% dplyr::select(sample_id, subject, age, time), by = "sample_id") %>%
  mutate(module_label = factor(display_label[module],
                               levels = display_label[non_grey]))

# Per-module × age deltas (for panel strip annotation)
delta_stats <- traj %>%
  group_by(module_label, age, subject) %>%
  summarize(
    pre  = if (any(time == "Pre"))  ME[time == "Pre"][1]  else NA_real_,
    post = if (any(time == "Post")) ME[time == "Post"][1] else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(delta = post - pre) %>%
  group_by(module_label, age) %>%
  summarize(mean_delta = mean(delta, na.rm = TRUE),
            sd_delta   = sd(delta, na.rm = TRUE), .groups = "drop") %>%
  mutate(stat_label = sprintf("\u0394\u0305 = %+.2f", mean_delta))

p <- ggplot(traj, aes(x = time, y = ME, group = subject, color = age)) +
  geom_line(alpha = 0.35, linewidth = 0.35) +
  geom_point(size = 0.9, alpha = 0.8) +
  stat_summary(aes(group = age), fun = mean, geom = "line",
               linewidth = 1.1, alpha = 0.95) +
  stat_summary(aes(group = age), fun = mean, geom = "point",
               size = 2.4, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.25) +
  geom_text(data = delta_stats,
            aes(x = 1.5, y = Inf, label = stat_label, color = age),
            vjust = 1.3, hjust = 0.5, inherit.aes = FALSE,
            size = 2.3, fontface = "bold",
            show.legend = FALSE) +
  facet_grid(age ~ module_label, scales = "free_y", switch = "y") +
  scale_color_manual(values = AGE_COLORS, guide = "none") +
  labs(title    = "Module Eigengene Trajectories (Pre \u2192 Post)",
       subtitle = "Per-subject ME paths | bold line = group mean | diamonds = group means",
       x = NULL, y = "Module eigengene") +
  FIG_THEME +
  theme(strip.text.x = element_text(size = 6, face = "bold", lineheight = 0.85),
        strip.text.y = element_text(size = 7.5, face = "bold"),
        strip.placement = "outside",
        panel.spacing.x = unit(1, "mm"),
        axis.text.x  = element_text(size = 6))

W <- 320; H <- 150
ggsave(file.path(RPT, "SUPP_me_trajectories.png"), p,
       width = W, height = H, units = "mm", dpi = 300)
ggsave(file.path(RPT, "SUPP_me_trajectories.pdf"), p,
       width = W, height = H, units = "mm", device = pdf_device)

message(sprintf("  SUPP ME trajectories saved (%d modules x 2 ages)",
                length(non_grey)))
