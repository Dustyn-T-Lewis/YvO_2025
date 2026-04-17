# Supplementary: WGCNA — Module Preservation (Pre -> Post Training)
# 200 permutations. Zsummary > 10: strong; 2-10: moderate; < 2: none.
# Generates: c_data/05_panel_E_preservation.csv (consumed by Panel A)
# NOTE: plot output removed — preservation visualized via Panel A Zsummary sidebar

setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(readr)
library(WGCNA)

DAT <- "04_Figures/F06/c_data"

meta           <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
datExpr        <- readRDS(file.path(DAT, "datExpr.rds"))
module_colors  <- readRDS(file.path(DAT, "module_colors.rds"))

mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
mod_display_vec    <- setNames(mod_bio_labels_df$display_label, mod_bio_labels_df$module_color)
mod_display_label  <- function(color) {
  lbl <- mod_display_vec[color]
  lbl[is.na(lbl)] <- stringr::str_to_title(color[is.na(lbl)])
  lbl
}

message("Panel E: module preservation...")

pre_samp_f  <- meta %>% filter(time == "Pre") %>% pull(sample_id)
post_samp_f <- meta %>% filter(time == "Post") %>% pull(sample_id)

datExpr_pre  <- datExpr[pre_samp_f, ]
datExpr_post <- datExpr[post_samp_f, ]

multiExpr <- list(
  Pre  = list(data = datExpr_pre),
  Post = list(data = datExpr_post)
)
multiColor <- list(Pre = module_colors)

cat("Starting module preservation (200 permutations, this may take 10-30 min)...\n")
mp <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  testNetworks      = 2,
  nPermutations     = 200,
  randomSeed        = 42,
  quickCor          = 0,
  verbose           = 3
)

ref  <- 1
test <- 2
z_summary <- mp$preservation$Z[[ref]][[test]]
mod_sizes <- z_summary[, "moduleSize"]

# Extract Z-component breakdown
z_density <- if ("Zdensity.pres" %in% colnames(z_summary)) {
  z_summary[, "Zdensity.pres"]
} else {
  rep(NA_real_, nrow(z_summary))
}
z_connectivity <- if ("Zconnectivity.pres" %in% colnames(z_summary)) {
  z_summary[, "Zconnectivity.pres"]
} else {
  rep(NA_real_, nrow(z_summary))
}

pres_df <- tibble(
  module           = rownames(z_summary),
  Zsummary         = z_summary[, "Zsummary.pres"],
  Zdensity         = z_density,
  Zconnectivity    = z_connectivity,
  module_size      = mod_sizes
) %>%
  filter(module != "gold", module != "grey") %>%
  mutate(preservation = case_when(
    Zsummary > 10 ~ "Strong",
    Zsummary > 2  ~ "Moderate",
    TRUE          ~ "Not preserved"
  ))

pres_df <- pres_df %>%
  mutate(bio_label = mod_display_label(module)) %>%
  arrange(Zsummary) %>%
  mutate(bio_label = factor(bio_label, levels = bio_label))

write_csv(pres_df, file.path(DAT, "05_panel_E_preservation.csv"))

message("  Panel E preservation CSV saved (plot output removed — consumed by Panel A)")
