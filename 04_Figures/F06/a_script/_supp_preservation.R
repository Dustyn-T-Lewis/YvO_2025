# Sourced by 02_supp_panels.R — expects style.R already loaded.

pacman::p_load(dplyr, readr, WGCNA)

BASE <- "04_Figures/F06"
DAT <- file.path(BASE, "c_data")

meta           <- read_csv(file.path(DAT, "meta.csv"))
datExpr        <- readRDS(file.path(DAT, "datExpr.rds"))
module_colors  <- readRDS(file.path(DAT, "module_colors.rds"))

mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"))
mod_display_vec    <- setNames(mod_bio_labels_df$display_label, mod_bio_labels_df$module_color)
mod_display_label  <- function(color) {
  lbl <- mod_display_vec[color]
  lbl[is.na(lbl)] <- stringr::str_to_title(color[is.na(lbl)])
  lbl
}

message("Panel E: module preservation...")

pre_samp_f  <- meta |> filter(time == "Pre") |> pull(sample_id)
post_samp_f <- meta |> filter(time == "Post") |> pull(sample_id)

datExpr_pre  <- datExpr[pre_samp_f, ]
datExpr_post <- datExpr[post_samp_f, ]

multiExpr <- list(
  Pre  = list(data = datExpr_pre),
  Post = list(data = datExpr_post)
)
multiColor <- list(Pre = module_colors)

pres_cache <- file.path(DAT, "modulePreservation_cache.rds")
if (file.exists(pres_cache)) {
  message("Loading cached modulePreservation result (delete to re-run)...")
  mp <- readRDS(pres_cache)
} else {
  message("Running modulePreservation (200 permutations, ~10-30 min)...")
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
  saveRDS(mp, pres_cache)
}

ref  <- 1
test <- 2
z_summary <- mp$preservation$Z[[ref]][[test]]
mod_sizes <- z_summary[, "moduleSize"]

pres_df <- tibble(
  module        = rownames(z_summary),
  Zsummary      = z_summary[, "Zsummary.pres"],
  Zdensity      = z_summary[, "Zdensity.pres"],
  Zconnectivity = z_summary[, "Zconnectivity.pres"],
  module_size   = mod_sizes
) |>
  filter(module != "gold", module != "grey") |>
  mutate(preservation = case_when(
    Zsummary > 10 ~ "Strong",
    Zsummary > 2  ~ "Moderate",
    TRUE          ~ "Not preserved"
  ))

pres_df <- pres_df |>
  mutate(bio_label = mod_display_label(module)) |>
  arrange(Zsummary) |>
  mutate(bio_label = factor(bio_label, levels = bio_label))

write_csv(pres_df, file.path(DAT, "05_panel_E_preservation.csv"))

message("  Panel E preservation CSV saved")
