# F02 Supp Panel C: CV Scatter Triptych (Pre vs Post + Δ Young vs Δ Old)
# C1/C2: per-protein CV% Pre vs Post (Young / Old)
# C3:    ΔCV Young vs ΔCV Old
# Outputs: pC (patchwork), SUPP_panel_C_cv_scatter.{pdf,png}

# Assumes style.R sourced, packages loaded, norm_df/norm_meta/samp_names set by parent

PA_SUB <- 60; PA_W <- 178; PA_H <- 70
RPT_PNG <- here::here("04_Figures", "F02", "b_reports", "supp", "png", "panels")
RPT_PDF <- here::here("04_Figures", "F02", "b_reports", "supp", "pdf", "panels")
DAT_DIR <- here::here("04_Figures", "F02", "c_data")

ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# Use norm_meta from parent (replaces regex parsing of sample IDs)
meta <- norm_meta
pdf_device <- get_pdf_device()

# CV on linear scale (Brenes 2024)
lin_mat <- 2^as.matrix(norm_df[, samp_names])

compute_cv <- function(mat, idx) {
  sub <- mat[, idx, drop = FALSE]
  apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
}

scatter_list <- lapply(c("Young", "Old"), function(ag) {
  pre_idx  <- meta$sample_id[meta$age == ag & meta$time == "Pre"]
  post_idx <- meta$sample_id[meta$age == ag & meta$time == "Post"]
  cv_pre   <- compute_cv(lin_mat, pre_idx)
  cv_post  <- compute_cv(lin_mat, post_idx)
  tibble(gene = norm_df$gene, cv_pre = cv_pre, cv_post = cv_post, age = ag)
})
scatter_df <- bind_rows(scatter_list) |>
  filter(!is.na(cv_pre), !is.na(cv_post)) |>
  mutate(
    delta_cv = cv_post - cv_pre,
    max_cv   = pmax(cv_pre, cv_post),
    age      = factor(age, levels = c("Young", "Old"))
  )

max_cv_cap <- quantile(scatter_df$max_cv, 0.98, na.rm = TRUE)
scatter_df$max_cv_capped <- pmin(scatter_df$max_cv, max_cv_cap)

cv_cap <- quantile(abs(scatter_df$delta_cv), 0.98, na.rm = TRUE)
scatter_df$delta_cv_capped <- pmin(pmax(scatter_df$delta_cv, -cv_cap), cv_cap)

top_cv_labels <- bind_rows(
  scatter_df |> filter(age == "Young") |> slice_max(max_cv, n = 15, with_ties = FALSE),
  scatter_df |> filter(age == "Old")   |> slice_max(max_cv, n = 15, with_ties = FALSE)
)

n_young <- sum(scatter_df$age == "Young" & !is.na(scatter_df$cv_pre) & !is.na(scatter_df$cv_post))
n_old   <- sum(scatter_df$age == "Old"   & !is.na(scatter_df$cv_pre) & !is.na(scatter_df$cv_post))
r_young <- cor(scatter_df$cv_pre[scatter_df$age == "Young"],
               scatter_df$cv_post[scatter_df$age == "Young"], use = "complete.obs")
r_old   <- cor(scatter_df$cv_pre[scatter_df$age == "Old"],
               scatter_df$cv_post[scatter_df$age == "Old"], use = "complete.obs")
ci_young <- fisher_z_ci(r_young, n_young)
ci_old   <- fisher_z_ci(r_old,   n_old)

r_annotations <- tibble(
  age   = factor(c("Young", "Old"), levels = c("Young", "Old")),
  label = c(sprintf("r = %.2f [%.2f, %.2f]", r_young, ci_young["lo"], ci_young["hi"]),
            sprintf("r = %.2f [%.2f, %.2f]", r_old,   ci_old["lo"],   ci_old["hi"]))
)

delta_wide <- scatter_df |>
  select(gene, delta_cv, age) |>
  pivot_wider(names_from = age, values_from = delta_cv, names_prefix = "dcv_") |>
  filter(!is.na(dcv_Young), !is.na(dcv_Old)) |>
  mutate(
    dist_origin     = sqrt(dcv_Young^2 + dcv_Old^2),
    mean_dcv        = (dcv_Young + dcv_Old) / 2,
    mean_dcv_capped = pmin(pmax(mean_dcv, -cv_cap), cv_cap)
  )

top_delta <- delta_wide |>
  slice_max(dist_origin, n = 15, with_ties = FALSE)

n_delta  <- sum(!is.na(delta_wide$dcv_Young) & !is.na(delta_wide$dcv_Old))
r_delta  <- cor(delta_wide$dcv_Young, delta_wide$dcv_Old, use = "complete.obs")
ci_delta <- fisher_z_ci(r_delta, n_delta)

theme_B <- FIG_THEME +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = FIG_STRIP_SIZE,
                                    face = "bold"),
    plot.subtitle    = element_text(size = FIG_SUBTITLE_SIZE - 1.0,
                                    face = "bold.italic", color = "grey40")
  )

axis_max_cv <- 300

pC12 <- ggplot(scatter_df, aes(x = cv_pre, y = cv_post)) +
  facet_wrap(~age, nrow = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = max_cv_capped), alpha = 0.35, size = 0.9) +
  geom_label_repel(data = top_cv_labels,
                   aes(label = gene, fill = max_cv_capped),
                   color = "white", fontface = "bold",
                   size = scale_text(BASE_GENE, PA_SUB),
                   label.padding = unit(1, "pt"),
                   linewidth = 0.3, max.overlaps = 20,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, seed = 42, show.legend = FALSE) +
  geom_label(data = r_annotations, aes(label = label),
             x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             size = scale_text(BASE_STAT, PA_SUB),
             color = "grey30", fontface = "bold",
             fill = alpha("white", 0.85), linewidth = 0,
             label.padding = unit(2, "pt"),
             inherit.aes = FALSE) +
  scale_color_viridis_c(option = "inferno", direction = -1,
                        begin = 0.1, end = 0.85, name = "CV%",
                        guide = guide_colorbar(barwidth = unit(2, "mm"),
                                               barheight = unit(12, "mm"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       begin = 0.1, end = 0.85, name = "CV%", guide = "none") +
  coord_fixed(xlim = c(0, axis_max_cv), ylim = c(0, axis_max_cv)) +
  labs(title = "Per-Protein Variability (CV%)",
       subtitle = sprintf(paste0("Per-protein CV%% Pre vs Post — measurement reproducibility | ",
                                 "%s proteins (cycloess-normalized)\n",
                                 "Young r = %.2f, Old r = %.2f | ",
                                 "delta-CV concordance r = %.2f (no systematic training shift)"),
                          format(nrow(norm_df), big.mark = ","),
                          r_young, r_old, r_delta),
       x = expression(bold(CV * "%"[Pre])),
       y = expression(bold(CV * "%"[Post])),
       tag = "c") +
  theme_B +
  theme(plot.title    = element_text(hjust = 0, size = FIG_TITLE_SIZE,
                                     face = "bold",
                                     margin = margin(b = 0)),
        plot.subtitle = element_text(hjust = 0, size = FIG_SUBTITLE_SIZE - 1.0,
                                     face = "bold.italic", color = "grey40",
                                     margin = margin(t = 0, b = 1)),
        strip.text    = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                     margin = margin(t = 0, b = 1)),
        legend.position      = c(0.97, 0.02),
        legend.justification = c(1, 0),
        legend.background    = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title         = element_text(face = "bold", size = FIG_LEGEND_TITLE),
        legend.text          = element_text(size = FIG_LEGEND_TEXT),
        legend.key.size      = unit(3, "mm"),
        plot.margin          = margin(t = 0, r = 0, b = 0, l = 5.5))

# Add a dummy facet column so "Training Response" renders as a strip header
# (matching the Young/Old facet strips in pC12 — same vertical position).
delta_wide$.facet <- "Training Response"

pC3 <- ggplot(delta_wide, aes(x = dcv_Young, y = dcv_Old)) +
  facet_wrap(~ .facet) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = mean_dcv_capped), alpha = 0.4, size = 0.9) +
  geom_label_repel(data = top_delta, aes(label = gene, fill = mean_dcv_capped),
                   color = "white", fontface = "bold",
                   size = scale_text(BASE_GENE, PA_SUB),
                   label.padding = unit(1, "pt"), linewidth = 0.3,
                   max.overlaps = 25,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, seed = 44, show.legend = FALSE) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("r = %.2f [%.2f, %.2f]",
                           r_delta, ci_delta["lo"], ci_delta["hi"]),
           hjust = -0.05, vjust = 1.4,
           size = scale_text(BASE_STAT, PA_SUB),
           color = "grey30", fontface = "bold",
           fill = alpha("white", 0.85), linewidth = 0,
           label.padding = unit(2, "pt")) +
  scale_color_gradient2(low = HEATMAP_LO, mid = "grey95", high = HEATMAP_HI,
                        midpoint = 0, limits = c(-cv_cap, cv_cap),
                        name = expression(bold(Delta * "CV%")),
                        guide = guide_colorbar(barwidth = unit(2, "mm"),
                                               barheight = unit(12, "mm"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  scale_fill_gradient2(low = HEATMAP_LO, mid = "grey95", high = HEATMAP_HI,
                       midpoint = 0, limits = c(-cv_cap, cv_cap),
                       guide = "none") +
  coord_fixed(xlim = c(-75, 225), ylim = c(-100, 200)) +
  labs(x = expression(bold(Delta * "CV%"[Young])),
       y = expression(bold(Delta * "CV%"[Old])),
       title = NULL) +
  theme_B +
  theme(strip.text           = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                            margin = margin(t = 0, b = 1)),
        legend.position      = c(0.97, 0.02),
        legend.justification = c(1, 0),
        legend.background    = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title         = element_text(face = "bold", size = FIG_LEGEND_TITLE),
        legend.text          = element_text(size = FIG_LEGEND_TEXT),
        legend.key.size      = unit(3, "mm"),
        axis.title.y         = element_text(margin = margin(r = 0, l = 0)),
        plot.margin          = margin(t = 0, r = 5.5, b = 0, l = 0))

write.csv(scatter_df |> select(gene, cv_pre, cv_post, delta_cv, age),
          file.path(DAT_DIR, "audit_panel_C_cv_scatter.csv"), row.names = FALSE)
write.csv(delta_wide |> select(gene, dcv_Young, dcv_Old, mean_dcv, dist_origin),
          file.path(DAT_DIR, "audit_panel_C_cv_delta.csv"), row.names = FALSE)

# 2:1 width ratio so the 2-facet (pC12) and 1-facet (pC3) sub-plots end up
# with equal-area panels. patchwork (vs cowplot::plot_grid) keeps pC as a
# proper ggplot/patchwork object so the supp stitch can override the panel
# tag via labs(tag = ...) without producing a duplicate tag.
pC <- (pC12 | pC3) + plot_layout(widths = c(2, 1))

ggsave(file.path(RPT_PNG, "SUPP_panel_C_cv_scatter.png"), pC,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_C_cv_scatter.pdf"), pC,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)

# --- Export for composite ---
pSA_title    <- "Per-Protein Variability (CV%)"
pSA_subtitle <- sprintf(
  "CV%% Pre vs Post | %s proteins | Y r = %.2f, O r = %.2f, \u0394CV r = %.2f",
  format(nrow(norm_df), big.mark = ","),
  r_young, r_old, r_delta)
pSA_legend   <- NULL
pC12 <- strip_for_composite(pC12)
pC3  <- strip_for_composite(pC3)

message("F02 Supp Panel C done")
