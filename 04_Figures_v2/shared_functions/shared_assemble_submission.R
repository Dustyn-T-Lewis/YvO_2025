library(here)
library(dplyr)
library(purrr)

manifest <- tribble(
  ~fig,  ~main_stem,            ~supp_stem,               ~fig_label,  ~supp_label, ~table_label,
  "F01", "MAIN_F01_composite",  "SUPP_F01_composite",     "Figure_1",  "S3_Figure", "S5_Table",
  "F02", "MAIN_F02_composite",  "SUPP_F02_composite",     "Figure_2",  "S4_Figure", "S6_Table",
  "F03", "MAIN_F03_composite",  "SUPP_F03_composite",     "Figure_3",  "S5_Figure", "S7_Table",
  "F04", "MAIN_F04_composite",  "SUPP_F04_diagnostics",   "Figure_4",  "S6_Figure", "S8_Table",
  "F05", "MAIN_F05_composite",  "SUPP_F05_diagnostics",   "Figure_5",  "S7_Figure", "S9_Table",
  "F06", "MAIN_F06_composite",  "SUPP_F06_composite",     "Figure_6",  "S8_Figure", "S10_Table",
  "F07", "MAIN_F07_composite",  "SUPP_F07_composite",     "Figure_7",  "S9_Figure", "S11_Table"
)

copy_one <- function(src, dest) {
  if (file.exists(src)) {
    file.copy(src, dest, overwrite = TRUE)
  } else {
    message("missing: ", src)
  }
}

supp_src <- function(dir, stem, ext) {
  exact <- file.path(dir, paste0(stem, ".", ext))
  if (file.exists(exact)) exact else file.path(dir, paste0(stem, "_main.", ext))
}

assemble_into <- function(pkg_dir) {
  fig_dir <- file.path(pkg_dir, "Main_Figures")
  sup_dir <- file.path(pkg_dir, "Supplementary")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

  pwalk(manifest, function(fig, main_stem, supp_stem, fig_label, supp_label, table_label, ...) {
    base <- here("04_Figures_v2", fig)
    for (ext in c("pdf", "png")) {
      copy_one(
        file.path(base, "b_reports", "main", ext, paste0(main_stem, ".", ext)),
        file.path(fig_dir, paste0(fig_label, ".", ext))
      )
      copy_one(
        supp_src(file.path(base, "b_reports", "supp", ext), supp_stem, ext),
        file.path(sup_dir, paste0(supp_label, ".", ext))
      )
    }
    copy_one(
      file.path(base, "c_data", paste0(fig, "_supplementary.xlsx")),
      file.path(sup_dir, paste0(table_label, ".xlsx"))
    )
  })

  copy_one(
    here("04_Figures_v2", "F00", "c_data", "F00_supplementary.xlsx"),
    file.path(sup_dir, "S4_Table.xlsx")
  )
}

walk(
  c("Submitted", "Rewrite"),
  ~ assemble_into(here("04_Figures", "Submission", .x))
)

message("submission figures assembled")
