#!/usr/bin/env Rscript
# Copy current pipeline outputs to the Box manuscript delivery folders.
# Refreshes everything the per-figure stitchers ship, without re-rendering.
# Set BOX to your manuscript root, then: Rscript sync_box.R

box <- Sys.getenv(
  "YVO_BOX",
  "~/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
)
box <- path.expand(box)

main_dir <- file.path(box, "02_Figures")
fig_dir <- file.path(box, "03_Supplementary", "figures")
tbl_dir <- file.path(box, "03_Supplementary", "tables")
for (d in c(
  file.path(main_dir, c("pdf", "png")),
  file.path(fig_dir, c("pdf", "png")), tbl_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# One row per delivered artifact. src is relative to the repo root (sans ext).
deliveries <- tibble::tribble(
  ~kind,  ~src,                                                       ~dest,
  "main", "04_Figures/F01/b_reports/main/%s/MAIN_F01_composite",      "MAIN_F01_composite",
  "main", "04_Figures/F02/b_reports/main/%s/MAIN_F02_composite",      "MAIN_F02_composite",
  "main", "04_Figures/F03/b_reports/main/%s/MAIN_F03_composite",      "MAIN_F03_composite",
  "main", "04_Figures/F04/b_reports/main/%s/MAIN_F04_composite",      "MAIN_F04_composite",
  "main", "04_Figures/F05/b_reports/main/%s/MAIN_F05_composite",      "MAIN_F05_composite",
  "main", "04_Figures/F06/b_reports/main/%s/MAIN_F06_composite",      "MAIN_F06_composite",
  "main", "04_Figures/F07/b_reports/main/%s/MAIN_F07_composite",      "MAIN_F07_composite",
  "fig",  "04_Figures/F00/b_reports/supp/%s/SUPP_F00_normalization",  "S01_Figure_Normalization",
  "fig",  "04_Figures/F00/b_reports/supp/%s/SUPP_F00_imputation",     "S02_Figure_Imputation",
  "fig",  "04_Figures/F01/b_reports/supp/%s/SUPP_F01_composite",      "S03_Figure_F01",
  "fig",  "04_Figures/F02/b_reports/supp/%s/SUPP_F02_composite",      "S04_Figure_F02",
  "fig",  "04_Figures/F03/b_reports/supp/%s/SUPP_F03_composite",      "S05_Figure_F03",
  "fig",  "04_Figures/F04/b_reports/supp/%s/SUPP_F04_diagnostics",    "S06_Figure_F04",
  "fig",  "04_Figures/F05/b_reports/supp/%s/SUPP_F05_diagnostics",    "S07_Figure_F05",
  "fig",  "04_Figures/F06/b_reports/supp/%s/SUPP_F06_composite",      "S08_Figure_F06",
  "tbl",  "01_normalization/c_data/01_normalization",                 "S01_Table_Normalization",
  "tbl",  "02_imputation/c_data/02_imputation",                       "S02_Table_Imputation",
  "tbl",  "03_DEP/c_data/03_DEP_results",                             "S03_Table_DEP",
  "tbl",  "04_Figures/F00/c_data/F00_supplementary",                  "S04_Table_F00",
  "tbl",  "04_Figures/F01/c_data/F01_supplementary",                  "S05_Table_F01",
  "tbl",  "04_Figures/F02/c_data/F02_supplementary",                  "S06_Table_F02",
  "tbl",  "04_Figures/F03/c_data/F03_supplementary",                  "S07_Table_F03",
  "tbl",  "04_Figures/F04/c_data/F04_supplementary",                  "S08_Table_F04",
  "tbl",  "04_Figures/F05/c_data/F05_supplementary",                  "S09_Table_F05",
  "tbl",  "04_Figures/F06/c_data/F06_supplementary",                  "S10_Table_F06",
  "tbl",  "04_Figures/F07/c_data/F07_supplementary",                  "S11_Table_F07"
)

copy_one <- function(src, dest_path) {
  if (!file.exists(src)) {
    warning("missing source: ", src, call. = FALSE)
    return(FALSE)
  }
  file.copy(src, dest_path, overwrite = TRUE)
}

ok <- 0L
for (i in seq_len(nrow(deliveries))) {
  d <- deliveries[i, ]
  if (d$kind == "tbl") {
    ok <- ok + copy_one(
      paste0(d$src, ".xlsx"),
      file.path(tbl_dir, paste0(d$dest, ".xlsx"))
    )
    next
  }
  out <- if (d$kind == "main") main_dir else fig_dir
  for (ext in c("pdf", "png")) {
    ok <- ok + copy_one(
      paste0(sprintf(d$src, ext), ".", ext),
      file.path(out, ext, paste0(d$dest, ".", ext))
    )
  }
}

# F07 supplementary ships the _main pdf and the plain composite png as S09.
ok <- ok + copy_one(
  "04_Figures/F07/b_reports/supp/pdf/SUPP_F07_composite_main.pdf",
  file.path(fig_dir, "pdf", "S09_Figure_F07.pdf")
)
ok <- ok + copy_one(
  "04_Figures/F07/b_reports/supp/png/SUPP_F07_composite.png",
  file.path(fig_dir, "png", "S09_Figure_F07.png")
)

message(sprintf("synced %d files to %s", ok, box))
