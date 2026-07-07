# Check bp_slim against a pinned GO Slim Generic release. Run on demand after
# editing bp_slim or to refresh the release; delete goslim_generic_bp.txt to refetch.

source("04_Figures_v2/shared/go_slim_terms.R")

release <- "2026-06-19"
cache <- "04_Figures_v2/shared/goslim_generic_bp.txt"
extras <- c("GO:0007586", "GO:0009100", "GO:0042180", "GO:0007126") # muscle additions + inert obsolete

if (!file.exists(cache)) {
  url <- sprintf("http://release.geneontology.org/%s/ontology/subsets/goslim_generic.obo", release)
  ids <- GSEABase::ids(GSEABase::getOBOCollection(url))
  ont <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, ids, "ONTOLOGY", "GOID"))
  bp <- sort(ont$GOID[!is.na(ont$ONTOLOGY) & ont$ONTOLOGY == "BP"])
  writeLines(c(
    sprintf("# GO Slim Generic biological-process subset, %d terms.", length(bp)),
    sprintf("# Source: %s", url),
    sprintf("# Pinned GO release %s; retrieved %s via GSEABase.", release, Sys.Date()),
    bp
  ), cache)
}

slim <- grep("^GO:", readLines(cache), value = TRUE)
unexplained <- setdiff(bp_slim, union(slim, extras))
if (length(unexplained)) {
  stop(
    "bp_slim terms absent from goslim_generic and undocumented: ",
    paste(unexplained, collapse = ", ")
  )
}
message(sprintf(
  "bp_slim OK: %d = %d goslim_generic BP + %d curated extras (release %s).",
  length(bp_slim), length(intersect(bp_slim, slim)),
  length(intersect(bp_slim, extras)), release
))
