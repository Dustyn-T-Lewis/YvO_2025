# Curated GO Slim biological-process terms and their 15-category consolidation.
# Single source of truth for go_slim_categories.R (gene->category assignment for
# the F04/F05 Panel B bars) and pathway_utils.R (build_goslim_gene_sets, fgsea/ORA).
#
# Provenance: 58 of the 62 terms are the GO Consortium goslim_generic
# biological-process subset (pinned release 2026-06-19, frozen in
# goslim_generic_bp.txt). Four are curated additions for muscle: GO:0007586
# digestion, GO:0009100 glycoprotein metabolism, GO:0042180 ketone metabolism, and
# GO:0007126 (meiosis, obsolete and inert). Run verify_go_slim.R to re-check the list
# against the pinned release.

bp_slim <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135"
)

SLIM_CONSOLIDATED <- c(
  "GO:0003012" = "Muscle & Contractile",
  "GO:0003013" = "Circulatory System",
  "GO:0030198" = "ECM & Adhesion", "GO:0007155" = "ECM & Adhesion",
  "GO:0034330" = "ECM & Adhesion", "GO:0042060" = "ECM & Adhesion",
  "GO:0007010" = "Cytoskeleton & Motility", "GO:0048870" = "Cytoskeleton & Motility",
  "GO:0007018" = "Cytoskeleton & Motility", "GO:0007163" = "Cytoskeleton & Motility",
  "GO:0044782" = "Cytoskeleton & Motility",
  "GO:0002376" = "Immune & Inflammation", "GO:0006954" = "Immune & Inflammation",
  "GO:0098542" = "Immune & Inflammation",
  "GO:0006629" = "Lipid Metabolism", "GO:0042180" = "Lipid Metabolism",
  "GO:0005975" = "Carbohydrate & Energy Metabolism",
  "GO:0006091" = "Carbohydrate & Energy Metabolism",
  "GO:1901135" = "Carbohydrate & Energy Metabolism",
  "GO:0006520" = "Amino Acid & Cofactor Metabolism",
  "GO:0055086" = "Amino Acid & Cofactor Metabolism",
  "GO:0006766" = "Amino Acid & Cofactor Metabolism",
  "GO:0071941" = "Amino Acid & Cofactor Metabolism",
  "GO:0098754" = "Amino Acid & Cofactor Metabolism",
  "GO:0007586" = "Amino Acid & Cofactor Metabolism",
  "GO:0007005" = "Mitochondria & Energy", "GO:0007031" = "Mitochondria & Energy",
  "GO:0006457" = "Protein Homeostasis", "GO:0030163" = "Protein Homeostasis",
  "GO:0006914" = "Protein Homeostasis", "GO:0051604" = "Protein Homeostasis",
  "GO:0065003" = "Protein Homeostasis", "GO:0009100" = "Protein Homeostasis",
  "GO:0055085" = "Transport", "GO:0016192" = "Transport",
  "GO:0006886" = "Transport", "GO:0006913" = "Transport",
  "GO:0072659" = "Transport", "GO:0061024" = "Transport",
  "GO:0002181" = "Translation & Ribosome", "GO:0042254" = "Translation & Ribosome",
  "GO:0006399" = "Translation & Ribosome",
  "GO:0006351" = "Transcription & Chromatin", "GO:0006355" = "Transcription & Chromatin",
  "GO:0016071" = "Transcription & Chromatin", "GO:0006325" = "Transcription & Chromatin",
  "GO:0006281" = "DNA & Cell Cycle", "GO:0006260" = "DNA & Cell Cycle",
  "GO:0006310" = "DNA & Cell Cycle", "GO:0032200" = "DNA & Cell Cycle",
  "GO:0000278" = "DNA & Cell Cycle", "GO:0140014" = "DNA & Cell Cycle",
  "GO:0007059" = "DNA & Cell Cycle", "GO:0000910" = "DNA & Cell Cycle",
  "GO:0007126" = "DNA & Cell Cycle",
  "GO:0048856" = "Development", "GO:0030154" = "Development",
  "GO:0012501" = "Development", "GO:0003014" = "Development",
  "GO:0003016" = "Development"
)

CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"              = "#E57373",
  "Cytoskeleton & Motility"           = "#FFB74D",
  "ECM & Adhesion"                    = "#FFF176",
  "Lipid Metabolism"                  = "#AED581",
  "Carbohydrate & Energy Metabolism"  = "#81C784",
  "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
  "Mitochondria & Energy"             = "#4DB6AC",
  "Protein Homeostasis"               = "#4FC3F7",
  "Transport"                         = "#7986CB",
  "Translation & Ribosome"            = "#BA68C8",
  "Transcription & Chromatin"         = "#AB47BC",
  "Immune & Inflammation"             = "#A1887F",
  "DNA & Cell Cycle"                  = "#90A4AE",
  "Circulatory System"                = "#CE93D8",
  "Development"                       = "#B0BEC5",
  "Other"                             = "#D0D0D0"
)
