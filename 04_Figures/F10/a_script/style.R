# F7 — Style
source("04_Figures/shared/style.R")

# --- F7 figure-level constants ---
F7_REF_W <- 200

F7_STRIP  <- scale_text(BASE_STAT, F7_REF_W)
F7_ANNOT  <- scale_text(BASE_GENE, F7_REF_W)
F7_TAG    <- F7_STRIP * 2
