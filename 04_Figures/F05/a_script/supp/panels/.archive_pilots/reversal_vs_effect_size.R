# F05 Supplementary: Reversal Probability vs Aging Effect Size
# Logistic regression + binned bars showing reversal rate by |logFC_Aging|.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)

RPT <- "04_Figures/F05/b_reports/supp/panels"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Load and prepare --------------------------------------------------------
res <- read.csv("03_DEP/c_data/03_combined_results.csv")

df <- res %>%
  select(gene, logFC_Aging, logFC_Training_Old, t_Aging) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  mutate(
    reversed       = as.integer(sign(logFC_Aging) != sign(logFC_Training_Old)),
    abs_logFC      = abs(logFC_Aging),
    abs_t          = abs(t_Aging),
    quintile_logFC = ntile(abs_logFC, 5),
    quintile_t     = ntile(abs_t, 5)
  )

# --- Binned reversal rates ---------------------------------------------------
bin_stats <- df %>%
  group_by(quintile_logFC) %>%
  summarise(
    n       = n(),
    n_rev   = sum(reversed),
    frac    = n_rev / n,
    mid_x   = median(abs_logFC),
    ci_lo   = binom.test(n_rev, n)$conf.int[1],
    ci_hi   = binom.test(n_rev, n)$conf.int[2],
    .groups = "drop"
  )

# --- Logistic regressions ----------------------------------------------------
fit_logfc <- glm(reversed ~ abs_logFC, data = df, family = binomial)
fit_t     <- glm(reversed ~ abs_t,     data = df, family = binomial)

or_logfc <- exp(coef(fit_logfc)["abs_logFC"])
ci_logfc <- exp(confint.default(fit_logfc)["abs_logFC", ])
p_logfc  <- summary(fit_logfc)$coefficients["abs_logFC", "Pr(>|z|)"]

or_t <- exp(coef(fit_t)["abs_t"])
ci_t <- exp(confint.default(fit_t)["abs_t", ])
p_t  <- summary(fit_t)$coefficients["abs_t", "Pr(>|z|)"]

# --- Logistic curve prediction -----------------------------------------------
pred_grid <- tibble(abs_logFC = seq(min(df$abs_logFC), max(df$abs_logFC),
                                    length.out = 200))
pred <- predict(fit_logfc, newdata = pred_grid, type = "link", se.fit = TRUE)
pred_grid <- pred_grid %>%
  mutate(
    fit  = plogis(pred$fit),
    lo   = plogis(pred$fit - 1.96 * pred$se.fit),
    hi   = plogis(pred$fit + 1.96 * pred$se.fit)
  )

# --- Save data ---------------------------------------------------------------
out <- bind_rows(
  bin_stats %>% mutate(type = "bin"),
  tibble(type = "logistic_logFC",
         or = or_logfc, or_lo = ci_logfc[1], or_hi = ci_logfc[2],
         p_value = p_logfc),
  tibble(type = "logistic_t",
         or = or_t, or_lo = ci_t[1], or_hi = ci_t[2],
         p_value = p_t)
)
write.csv(out, file.path(DAT, "supp_reversal_vs_effect.csv"),
          row.names = FALSE)

# --- Visualization -----------------------------------------------------------
p <- ggplot() +
  # Logistic curve + ribbon
  geom_ribbon(data = pred_grid, aes(x = abs_logFC, ymin = lo, ymax = hi),
              fill = "#7B5EA7", alpha = 0.15) +
  geom_line(data = pred_grid, aes(x = abs_logFC, y = fit),
            color = "#7B5EA7", linewidth = 0.8) +
  # Binned bars
  geom_point(data = bin_stats, aes(x = mid_x, y = frac),
             size = 3, color = "#2E7D32") +
  geom_errorbar(data = bin_stats, aes(x = mid_x, ymin = ci_lo, ymax = ci_hi),
                width = 0.02, linewidth = 0.5, color = "#2E7D32") +
  geom_text(data = bin_stats,
            aes(x = mid_x, y = ci_hi + 0.02,
                label = sprintf("%d/%d", n_rev, n)),
            size = 2.5, color = "grey30") +
  # Reference line
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  # Annotation
  annotate("label", x = max(df$abs_logFC) * 0.95, y = 0.15,
           hjust = 1, size = 2.8, fill = "white", alpha = 0.85,
           label.size = 0.3,
           label = sprintf("logFC model: OR = %.2f [%.2f, %.2f], %s\nt-stat model: OR = %.2f [%.2f, %.2f], %s",
                           or_logfc, ci_logfc[1], ci_logfc[2], fmt_p(p_logfc),
                           or_t, ci_t[1], ci_t[2], fmt_p(p_t))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title    = "Reversal Probability vs Aging Effect Size",
    subtitle = "Binned reversal rates (green) with logistic fit (purple)",
    x = "|logFC(Aging)|",
    y = "Reversal probability"
  ) +
  FIG_THEME

ggsave(file.path(RPT, "SUPP_reversal_vs_effect.png"), p,
       width = 200, height = 150, units = "mm", dpi = 300)

cat("Reversal vs effect size plot saved.\n")
cat(sprintf("  logFC model: OR = %.3f [%.3f, %.3f], %s\n",
            or_logfc, ci_logfc[1], ci_logfc[2], fmt_p(p_logfc)))
cat(sprintf("  t-stat model: OR = %.3f [%.3f, %.3f], %s\n",
            or_t, ci_t[1], ci_t[2], fmt_p(p_t)))
